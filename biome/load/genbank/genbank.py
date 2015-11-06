from Bio import SeqIO, SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from py2neo import neo4j, node, rel, cypher
from time import ctime, time
import warnings
import logging
import os


def parse2read(gb_file):
    multi_file = SeqIO.parse(gb_file, 'genbank')
    return [rec for rec in multi_file]

def node2link(gb_node):
    return str(gb_node).split(' ')[0].split('(')[1]

def num2strand(num):
    if num == 1:
        strand = 'forward'
    elif num == -1:
        strand = 'reverse'
    else:
        strand = 'unknown'
    return strand

class BioGraphConnection():
    def __init__(self, db_link='http://localhost:8484/db/data/'):
        if not isinstance(db_link, str):
            raise TypeError ('Argument must be string.')

        self.db_link = db_link
        self.data_base = neo4j.GraphDatabaseService(self.db_link)

class GenBank():
    non_gene_list = ['mobile_element', 'repeat_region', 'rep_origin', 'STS']
    gene_product_list = ['CDS', 'tRNA', 'rRNA', 'tmRNA', 'ncRNA', 'mRNA']
    def __init__(self, gb_file, db_connection, logger_level = logging.INFO):
        if not isinstance(gb_file, str):
            raise TypeError('gb_file must be a string with filename.')
        if not isinstance(db_connection, BioGraphConnection):
            raise TypeError('db_connection must be an instance of the class BioGraphConnection.')
        if not os.path.isfile(gb_file):
            raise ValueError('There is no %s in current directory.' % gb_file)
        logging.basicConfig(filename = 'BiomeDB.log', level = logger_level,
                            format = '%(asctime)s %(message)s - %(module)s.%(funcName)s',
                            datefmt='%H:%M:%S-%d.%m.%y')
        self._logger1 = logging.getLogger('GenBank.log')
        # self._logger1 = logging.getLogger(__name__)
        self._logger1.info('GenBank object was created. For file ' + gb_file)
        self.gb_file = gb_file
        self.rec, self.seq_type = self._read_gb_file()
        self.db_connection = db_connection
        self.organism_list = [None, None, None]
        self.ccp_list = [None, None, None]
        self.external_sources = self.get_external_db()

    def get_external_db(self):
        db_dict = {}
        dbs = list(self.db_connection.data_base.find('DB'))
        for db in dbs:
            db_dict[db.get_properties()['name']] = db
        self._logger1.info('The list of external data bases was obtained.')
        if not db_dict.has_key('GenBank'):
            genbank_node, = self.db_connection.data_base.create(node({'name':'GenBank'}))
            genbank_node.add_labels('DB')
            db_dict['GenBank'] = genbank_node
            self._logger1.info('The GenBank node was not found and was created.')
        return db_dict

    def upload(self):
        """

        """

        # Method searches pattern (org)<-[:PART_OF]-(ccp)-[:EVIDENCE]->(xref).
        # If it is not found it searches for organism, if organism is not found,
        # it is created. The ccp is being searched, if it is not found, it is created
        # with its xref.
        self.search_pattern_cypher()

        # Feature loop.
        self._logger.info('Start reading the gb-file features.')
        for i, feature in enumerate(self.rec.features):
            if feature.type == 'gene':
                self.make_gene_and_product(i)
            elif feature.type in self.non_gene_list or feature.type[:4] == 'misc':
                self.make_non_gene(feature)
            elif feature.type in self.gene_product_list:
                pass
            else:
                self._logger.info('Unknown element %s was skipped.' % feature.type)

    def _read_gb_file(self):
        rec = SeqIO.read(self.gb_file, 'genbank')
        if not rec:
            self._logger.error('gb-file is empty.')
            raise ValueError('gb-file is empty.')
        seq_file = open(self.gb_file, 'r')
        first_line = seq_file.readline()
        seq_file.close()
        if 'circular' in first_line:
            seq_type = 'circular'
        else:
            seq_type = 'linear'
        self._logger.info('The header of gb-file was read successfully.')
        return rec, seq_type

    def search_pattern_cypher(self):
        self._logger.info('Searching for start pattern:')
        taxon = self.rec.features[0].qualifiers['db_xref'][0].split(':')[1]
        session = cypher.Session(self.db_connection.db_link)
        transaction = session.create_transaction()
        query = 'MATCH (org:Organism)<-[:PART_OF]-(ccp)-[:EVIDENCE]->(xref:XRef) ' \
                'WHERE org.name="%s" AND ccp.length=%d AND xref.id="%s" RETURN org, ccp, xref' \
                % (self.rec.annotations['source'],
                   int(self.rec.features[0].location.end),
                   self.rec.annotations['accessions'][0])
        transaction.append(query)
        transaction_res = list(transaction.commit())[0]
        self._logger.info(query)
        if transaction_res:
            self._logger.info('Pattern was found.')
            transaction_res = transaction_res[0]
            self.organism_list = [transaction_res[0], self.rec.annotations['organism'], node2link(transaction_res[0])]
            self.update_source(transaction_res[0], transaction_res[0].get_properties())
            self.ccp_list = [transaction_res[1], self.rec.description, node2link(transaction_res[1])]
            self.update_source(transaction_res[1], transaction_res[1].get_properties())
        else:
            self._logger1.info('Pattern was not found.')
            # Create or find organism
            self.create_or_update_organism()

            # Create or find nodes for chromosome, contig or plasmid
            self.create_or_update_ccp()

    def create_or_update_organism(self):
        # searching organism
        organism_name = self.rec.annotations['organism']
        self._logger1.info('Searching for organism %s' % organism_name)
        search_organism = list(self.db_connection.data_base.find('Organism', 'name', organism_name))
        if not search_organism:
            self._logger1.info('Organism was not found, creating organism node.')
            # creating organism
            current_organism, = self.db_connection.data_base.create(
                node({'name': organism_name,
                      'source': ['GenBank']}))
            current_organism.add_labels('Organism')
        elif len(search_organism) > 1:
            self._logger1.warning('There are duplicates in the DB.')
            raise Exception('There are duplicates in the DB.')
        else:
            self._logger1.info('Organism was found.')
            current_organism = search_organism[0]
            self.update_source(current_organism, current_organism.get_properties())
        self.organism_list = [current_organism, organism_name, node2link(current_organism)]

    def create_or_update_ccp(self):
        description = self.rec.description
        if 'complete genome' in description or 'complete sequence' in description\
                and not 'lasmid' in description:
            ccp_label = 'Chromosome'
        elif 'ontig' in description:
            ccp_label = 'Contig'
        elif 'lasmid' in description:
            ccp_label = 'Plasmid'
        else:
            self._logger1.error('Unknown genome element')
            warnings.warn('Unknown genome element')
            ccp_label = 'Contig'
        ccp_length = len(self.rec)
        ccp_name = description
        self._logger1.info('Searching for %s with name %s.' % (ccp_label, ccp_name))
        search_ccp = self.search_node(ccp_label, ['name', 'length'],
                                      [ccp_name, ccp_length])
        if not search_ccp:
            self._logger1.info('%s was not found.' % ccp_label)
            # Creating chromosome, contig or plasmid
            current_ccp, = self.db_connection.data_base.create(node({'name': ccp_name,
                                                                   'length': ccp_length,
                                                                   'type': self.seq_type,
                                                                   'source': ['GenBank']}))
            # Adding label
            current_ccp.add_labels(ccp_label, 'BioEntity', 'DNA')
            self._logger1.info('%s was created.' % ccp_label)
        else:
            current_ccp = search_ccp[0][0]
            self.update_source(current_ccp, current_ccp.get_properties())
            self._logger1.info('%s was found.' % ccp_label)
        self.ccp_list = [current_ccp, ccp_name, node2link(current_ccp)]
        self.db_connection.data_base.create(rel(self.ccp_list[0], 'PART_OF', self.organism_list[0]))

        # Check XRefs
        refseq = self.rec.annotations['accessions'][0]
        self._logger1.info('Searching for XRefs and %s pattern:' % ccp_label)
        session = cypher.Session(self.db_connection.db_link)
        transaction = session.create_transaction()
        query = 'START ccp=node(%s) MATCH (ccp)-[e:EVIDENCE]->(xref) WHERE xref.id="%s" RETURN e' \
                % (self.ccp_list[2], refseq)
        self._logger1.info(query)
        transaction.append(query)
        transaction_res = list(transaction.commit())[0]
        if not transaction_res:
            self._logger1.info('Pattern was not found. Searching for XRef node.')
            check_xref = self.search_node('XRef', ['id'], [refseq])
            if not check_xref:
                # Create XRefs
                self._logger1.info('XRef node was not found. Creating XRef node.')
                refseq_node, evidence_link, refseq_link = self.db_connection.data_base.\
                    create(node({'id': refseq}),
                           rel(current_ccp, 'EVIDENCE', 0),
                           rel(0, 'LINK_TO', self.external_sources['GenBank']))
                refseq_node.add_labels('XRef')
            else:
                self._logger1.info('XRef node was found. Creating LINK_TO relation.')
                # Check relation
                evidence, link_to = self.db_connection.data_base.\
                    create(rel(current_ccp, 'EVIDENCE', check_xref[0][0]),
                           rel(check_xref[0][0], 'LINK_TO', self.external_sources['GenBank']))
        else:
            self._logger1.info('Pattern was found.')


    def make_gene_and_product(self, i):
        self._logger1.info('Processing gene.')
        gene = self.rec.features[i]
        # Take info out of gene
        gene_dict = {'locus_tag':gene.qualifiers['locus_tag'][0]}
        try:
            next_feature = self.rec.features[i+1]
            if not next_feature.type in self.gene_product_list:
                next_feature = self.rec.features[i]
                self._logger1.info('Gene has no product.')
            else:
                gene_dict['product'] = next_feature.qualifiers['product'][0]
                self._logger1.info('Gene has product %s' % next_feature.type)
                # pass
        except:
            next_feature = self.rec.features[i]
            self._logger1.info('Gene product finding fails.')
        gene_dict['start'] = int(next_feature.location.start) + 1
        gene_dict['end'] = int(next_feature.location.end)
        gene_dict['strand'] = num2strand(next_feature.location.strand)

        # Setting the gene name
        try:
            gene_name = next_feature.qualifiers['gene'][0]
        except:
            gene_name = next_feature.qualifiers['locus_tag'][0]
        gene_dict['name'] = gene_name

        gene_node, check = self.create_or_update_gene(gene_dict)

        # Creating XRef
        try:
            xrefs = gene.qualifiers['db_xref']

            # Removing GIs from xref list
            for xref in xrefs:
                if 'GI' in xref:
                    xrefs.pop(xref)

            self.create_xref(xrefs=xrefs, feature_node=gene_node, check=check)
            self._logger1.info('Gene xrefs were found.')
        except:
            self._logger1.info('Gene has no xrefs.')

        if next_feature.type == 'CDS':
            try:
                self.create_or_update_cds(next_feature, gene_node)
            except ValueError:
                pass
        elif next_feature.type in self.gene_product_list[1:-1]:
            self.create_or_update_rna(next_feature, gene_node)
        else:
            pass


    def create_or_update_gene(self, gene_dict):
        # Searching the gene in the DB
        self._logger1.info('Searching gene %s' % gene_dict['name'])
        check_gene = self.search_gene_pattern(gene_dict['start'],
                                              gene_dict['end'],
                                              gene_dict['strand'])

        # If gene is not found
        if not check_gene:
            self._logger1.info('Gene was not found. Creating gene.')
            # Creating a new gene
            gene_dict['source'] = ['GenBank']
            gene_node, part_of_org, part_of_ccp = self.db_connection.data_base.create(
                    node(gene_dict),
                    rel(0, 'PART_OF', self.organism_list[0]),
                    rel(0, 'PART_OF', self.ccp_list[0]))

            # Creating labels for a new gene node
            gene_node.add_labels('Gene', 'Feature', 'DNA', 'BioEntity')

            # Creating term
            self.create_term(gene_dict['name'], gene_node)

            # Flag for create_xref method.
            updated = False
        else:
            self._logger1.info('Gene was found. Updating gene.')

            # Updating the gene
            gene_node = check_gene[0][0]
            gene_props = gene_node.get_properties()

            # Updating source
            self.update_source(gene_node, gene_props)

            # Checking locus_tag
            try:
                self._logger1.info('Checking locus_tag.')
                gene_locus_tag = gene_props['locus_tag']
                if gene_locus_tag != gene_dict['locus_tag']:
                    self._logger1.info('Locus_tag is different.')
                else:
                    self._logger1.info('Locus_tag matches.')
            except:
                gene_node.update_properties({'locus_tag': gene_dict['locus_tag']})
                self._logger1.info('Locus_tag was created.')

            # Checking product
            try:
                self._logger1.info('Checking product.')
                gene_product = gene_props['product']
                if gene_product != gene_dict['product']:
                    self._logger1.info('Product is different.')
                else:
                    self._logger1.info('Product matches.')
            except:
                try:
                    gene_node.update_properties({'product': gene_dict['product']})
                    self._logger1.info('Product was created.')
                except:
                    self._logger1.info('Product was not created.')

            # Flag for create_xref method.
            updated = True

            # Checking name of existing gene
            if gene_props['name'] != gene_dict['name']:
                self._logger1.info('Creating additional Term for gene.')
                self.create_term(gene_dict['name'], gene_node)
        return gene_node, updated

    def create_or_update_rna(self, feature, gene_node):
        # Taking info for the RNA
        self._logger1.info('Processing RNA.')
        try:
            rna_dict = {'name':feature.qualifiers['product'][0]}
        except:
            rna_dict = {'name':feature.qualifiers['locus_tag'][0]}

        try:
            rna_dict['comment'] = feature.qualifiers['note'][0]
        except:
            pass

        # Searching for the RNA
        self._logger1.info('Searching for RNA:')
        session = cypher.Session(self.db_connection.db_link)
        transaction = session.create_transaction()
        query = 'START g=node(%s) ' \
                'MATCH (g)-[:ENCODES]->(r:RNA) ' \
                'RETURN r' % node2link(gene_node)
        self._logger1.info(query)
        transaction.append(query)
        check_rna = transaction.commit()[0]
        if not check_rna:
            self._logger1.info('RNA was not found. Creating RNA.')
            rna_dict['source'] = ['GenBank']
            # Creating a RNA node
            rna_node, part_of_org, encodes = self.db_connection.data_base.create(
                    node(rna_dict),
                    rel(0, 'PART_OF', self.organism_list[0]),
                    rel(gene_node, 'ENCODES', 0))

            # Renaming ncRNA to sRNA
            if feature.type == 'ncRNA':
                feature.type = 'sRNA'

            # Adding labels to the node
            rna_node.add_labels(feature.type, 'RNA', 'BioEntity')

            # Creating term
            self.create_term(rna_dict['name'], rna_node)

            # Flag for create_xref method.
            updated = False
        else:
            # Updating the RNA node
            self._logger1.info('RNA was found. Updating RNA.')
            rna_node = check_rna[0][0]
            rna_props = rna_node.get_properties()

            # Updating source
            self.update_source(rna_node, rna_props)

            # Checking name of existing gene
            if rna_props['name'] != rna_dict['name']:
                self.create_term(rna_dict['name'], rna_node)

            # Flag for create_xref method.
            updated = True

        # Creating XRef
        try:
            self._logger1.info('Creating XRefs for RNA')
            xrefs = [db for db in feature.qualifiers['db_xref'] if db.split(':')[0] not in ['GeneID']]
            self.create_xref(xrefs=xrefs, feature_node=rna_node, check=updated),
        except:
            self._logger1.info('RNA has no xrefs.')

    def create_or_update_cds(self, feature, gene_node):
        # Taking info for the polypetide
        self._logger1.info('Processing CDS.')
        poly_dict = {'name':feature.qualifiers['product'][0]}
        try:
            seq = feature.qualifiers['translation'][0]
            self._logger1.info('CDS has translation in source.')
        except:
            dna_seq = self.rec.seq[feature.location.nofuzzy_start:feature.location.nofuzzy_end]
            try:
                seq = dna_seq.translate(feature.qualifiers['trans_table'][0])[:-1]
                self._logger1.info('CDS has no translation in source. Translated manually.')
            except:
                self._logger1.info('Gene was not translated. start: %d, end: %d, seq: %s' % (feature.location.nofuzzy_start, feature.location.nofuzzy_end, dna_seq))
                raise ValueError
        poly_dict['seq'] = seq
        try:
            poly_dict['comment'] = feature.qualifiers['note'][0]
            self._logger1.info('CDS has comment.')
        except:
            self._logger1.info('CDS has no comment.')

        # Searching for the polypeptide
        self._logger1.info('Searching for CDS:')
        session = cypher.Session(self.db_connection.db_link)
        transaction = session.create_transaction()
        query = 'START g=node(%s) ' \
                'MATCH (g)-[:ENCODES]->(p:Polypeptide) ' \
                'WHERE p.seq="%s" ' \
                'RETURN p' % (node2link(gene_node), seq)
        self._logger1.info(query)
        transaction.append(query)
        check_poly = transaction.commit()[0]
        if not check_poly:
            self._logger1.info('CDS was not found. Creating CDS.')
            poly_dict['source'] = ['GenBank']
            # %%%%%%%%%%%
            if 'hypothetical protein' in poly_dict['name']:
                poly_dict['name'] = 'hypothetical protein %s' %feature.qualifiers['locus_tag'][0]
            # %%%%%%%%%%%
            # Creating a polypeptide node
            poly_node, part_of_org, encodes = self.db_connection.data_base.create(
                    node(poly_dict),
                    rel(0, 'PART_OF', self.organism_list[0]),
                    rel(gene_node, 'ENCODES', 0))

            # Adding labels to the node
            poly_node.add_labels('Polypeptide', 'Peptide', 'BioEntity')

            # Creating term
            self.create_term(poly_dict['name'], poly_node)

            # Flag for create_xref method.
            updated = False
        else:
            # Updating the polypeptide node
            self._logger1.info('CDS was found. Updating CDS.')
            poly_node = check_poly[0][0]
            poly_props = poly_node.get_properties()

            # Updating source
            self.update_source(poly_node, poly_props)

            # Checking name of existing gene
            if poly_props['name'] != poly_dict['name']:
                self.create_term(poly_dict['name'], poly_node)

            # Flag for create_xref method.
            updated = True

        # Creating XRef
        try:
            self._logger1.info('Creating XRefs for CDS')
            xrefs = [db for db in feature.qualifiers['db_xref'] if db.split(':')[0] not in ['GeneID']]
            self.create_xref(xrefs=xrefs, feature_node=poly_node, check=updated),
        except:
            self._logger1.info('CDS has no xrefs.')

    def make_non_gene(self, feature):
        self._logger1.info('Processing non-gene feature.')
        feature_dict ={}
        feature_dict['start'] = int(feature.location.start) + 1
        feature_dict['end'] = int(feature.location.end)
        feature_dict['strand'] = num2strand(feature.location.strand)
        check_feature = self.search_non_gene_pattern(feature_dict['start'],
                                              feature_dict['end'],
                                              feature_dict['strand'])

        self._logger1.info('Searching for the  feature.')

        try:
            feature_dict['comment'] = feature.qualifiers['note'][0]
            self._logger1.info('Feature has comment.')
        except:
            self._logger1.info('Feature has no comment.')
        try:
            feature_dict['type'] = feature.qualifiers['mobile_element_type']
        except:
            pass
        if not check_feature:
            self._logger1.info('Feature was not found. Creating feature.')
            feature_dict['source'] = ['GenBank']
            # Creating a feature node
            feature_node, part_of_org = self.db_connection.data_base.create(
                    node(feature_dict),
                    rel(0, 'PART_OF', self.ccp_list[0]))
            if feature.type in self.non_gene_list[:-1]:
                feature.type = feature.type.capitalize()
            feature_node.add_labels(feature.type, 'Feature', 'DNA')

            # Flag for create_xref method.
            updated = False
        else:
            self._logger1.info('Feature was found. Updating feature.')
            # Updating the feature
            feature_node = check_feature[0][0]
            feature_props = feature_node.get_properties()

            # Updating source
            self.update_source(feature_node, feature_props)

            # Flag for create_xref method.
            updated = True

            # Checking comment of the existing feature.
            if feature_dict.has_key('comment'):
                feature_node.update_properties({'comment': feature_dict['comment']})
        # Creating xrefs
        try:
            self._logger1.info('Creating XRefs for feature.')
            xrefs = feature.qualifiers['db_xref']
            self.create_xref(xrefs=xrefs, feature_node=feature_node, check=updated),
        except:
            self._logger1.info('Feature has no xrefs.')

    def update_source(self, node, props):
        try:
            source = props['source']
            if not isinstance(source, list):
                source = [source]
            if not 'GenBank' in source:
                source.append('GenBank')
                node.update_properties({'source': sorted(source)})
        except:
            node.update_properties({'source': ['GenBank']})

    # def search_gene_pattern(self, start, end, strand):
    #     if not isinstance(start, int) or not isinstance(end, int):
    #         raise ValueError('start and end must be string.')
    #     if not isinstance(strand, str):
    #         raise ValueError('strand must be string.')
    #     if not strand in ["forward", "reverse", "unknown"]:
    #         raise ValueError('strand must be "forward", "reverse" or "unknown".')
    #
    #     self._logger1.info('Searching gene pattern:')
    #     session = cypher.Session(self.db_connection.db_link)
    #     transaction = session.create_transaction()
    #     query = 'START org=node(%s), ccp=node(%s) ' \
    #             'MATCH (ccp)<-[:PART_OF]-(g:Gene)-[:PART_OF]->(org) ' \
    #             'WHERE g.start=%d AND g.end=%d AND g.strand="%s" ' \
    #             'RETURN g' \
    #             % (self.organism_list[2], self.ccp_list[2], start, end, strand)
    #     self._logger1.info(query)
    #     transaction.append(query)
    #     try:
    #         transaction_res = list(transaction.commit())[0]
    #         self._logger1.info('Transaction succeeded.')
    #         return transaction_res
    #     except:
    #         self._logger1.error('Transaction failed.')

    def search_gene_pattern(self, start, end, strand):
        if not isinstance(start, int) or not isinstance(end, int):
            raise ValueError('start and end must be string.')
        if not isinstance(strand, str):
            raise ValueError('strand must be string.')
        if not strand in ["forward", "reverse", "unknown"]:
            raise ValueError('strand must be "forward", "reverse" or "unknown".')

        self._logger1.info('Searching gene pattern "ORGANISM--GENE--CCP":')
        session = cypher.Session(self.db_connection.db_link)
        transaction = session.create_transaction()
        try:
            query = 'START org=node(%s), ccp=node(%s) ' \
                    'MATCH (ccp)<-[:PART_OF]-(g:Gene)-[:PART_OF]->(org) ' \
                    'WHERE g.start=%d AND g.end=%d AND g.strand="%s" ' \
                    'RETURN g' \
                    % (self.organism_list[2], self.ccp_list[2], start, end, strand)
            self._logger1.info(query)
            transaction.append(query)
            transaction_res = list(transaction.commit())[0]
            self._logger1.info('Transaction succeeded.')
            if not transaction_res:
                self._logger1.error('Transaction failed.')
                try:
                    self._logger1.info('Searching gene pattern "GENE--POLYPEPTIDE--ORGANISM":')
                    session = cypher.Session(self.db_connection.db_link)
                    transaction = session.create_transaction()
                    query = 'START org=node(%s) ' \
                            'MATCH (g:Gene)-[:ENCODES]->(p:Polypeptide)-[:PART_OF]->(org) ' \
                            'WHERE g.start=%d AND g.end=%d AND g.strand="%s" ' \
                            'RETURN g, p' \
                            % (self.organism_list[2], start, end, strand)
                    self._logger1.info(query)
                    transaction.append(query)
                    transaction_res = list(transaction.commit())[0]
                    self._logger1.info('Transaction succeeded.')
                    print 'Make edge'
                    self.db_connection.data_base.create(rel(transaction_res[0][0], 'PART_OF', self.organism_list[0]),
                                                        rel(transaction_res[0][0], 'PART_OF', self.ccp_list[0]))
                    print 'Edge is made'
                except:
                    pass
            return transaction_res
        except:
            self._logger1.error('Transaction failed.')

    def search_non_gene_pattern(self, start, end, strand):
        if not isinstance(start, int) or not isinstance(end, int):
            raise ValueError('start and end must be string.')
        if not isinstance(strand, str):
            raise ValueError('strand must be string.')
        if not strand in ["forward", "reverse", "unknown"]:
            raise ValueError('strand must be "forward", "reverse" or "unknown".')

        self._logger1.info('Searching non-gene pattern:')
        session = cypher.Session(self.db_connection.db_link)
        transaction = session.create_transaction()
        query = 'START ccp=node(%s) ' \
                'MATCH (ccp)<-[:PART_OF]-(ng) ' \
                'WHERE ng.start=%d AND ng.end=%d AND ng.strand="%s" ' \
                'RETURN ng' \
                % (self.ccp_list[2], start, end, strand)
        self._logger1.info(query)
        transaction.append(query)
        try:
            transaction_res = list(transaction.commit())[0]
            self._logger1.info('Transaction succeeded.')
            return transaction_res
        except:
            self._logger1.error('Transaction failed.')

    def create_term(self, name, bioentity):
        if not isinstance(name, str):
            raise ValueError('name must be string.')
        if not isinstance(bioentity, neo4j.Node):
            raise ValueError('bioentity must be py2neo.neo4j.Node.')

        # Searching term
        self._logger1.info('Searching term pattern:')
        check_term = list(self.db_connection.data_base.find('Term', 'text', name))
        if not check_term:
            term, has_name = self.db_connection.data_base.create(node({'text': name}),
                                                                 rel(bioentity, 'HAS_NAME', 0))
            term.add_labels('Term')
            self._logger1.info('Term and "HAS_NAME" relation were created.')
        else:
            if len(check_term) == 1:
                session = cypher.Session(self.db_connection.db_link)
                transaction = session.create_transaction()
                query = 'START entity=node(%s), term=node(%s) ' \
                        'MATCH (entity)-[r:HAS_NAME]->(term) ' \
                        'RETURN r' \
                        % (node2link(check_term[0]), node2link(bioentity))
                self._logger1.info(query)
                transaction.append(query)
                transaction_res = list(transaction.commit())[0]
                if not transaction_res:
                    has_name, = self.db_connection.data_base.create(rel(bioentity, 'HAS_NAME', check_term[0]))
                    self._logger1.info('Term was found and "HAS_NAME" relation was created.')
                else:
                    if len(transaction_res[0]) == 1:
                        pass
                    else:
                        self._logger1.info('There is a duplicate "HAS_NAME" relation with Term: %s' % name)
            else:
                has_name, = self.db_connection.data_base.create(rel(bioentity, 'HAS_NAME', check_term[0]))
                self._logger1.info('There is a duplicate Term: %s' % name)

        # self.db_connection.data_base.create(rel(bioentity, 'EVIDENCE', xref))

    def create_xref(self, xrefs, feature_node, check=False):
        """
        Method creates XRef to a feature with relation EVIDENCE.
        :param xrefs: list
        :param feature_node: py2neo.neo4j.Node
        :return: None
        """
        if not isinstance(xrefs, list):
            raise ValueError('refs must be a list.')
        if not isinstance(feature_node, neo4j.Node):
            raise ValueError('feature must be py2neo.neo4j.Node.')

        self._logger1.info('Creating XRefs.')
        id_list = []

        # Check if feature is already has an XRef
        if check == True:
            self._logger1.info('Searching for XRef pattern:')
            session = cypher.Session(self.db_connection.db_link)
            transaction = session.create_transaction()
            query = 'START feature=node(%s)' \
                    ' MATCH (feature)-[:EVIDENCE]->(xref:XRef)-[:LINK_TO]->(db)' \
                    ' RETURN xref.id, db.name, db' % node2link(feature_node)
            self._logger1.info(query)
            transaction.append(query)
            transaction_res = list(transaction.commit())[0]
            for rec in transaction_res:
                id_list.append(rec[0])
                if not self.external_sources.has_key(rec[1]):
                    self.external_sources[rec[1]] = rec[2]
            self._logger1.info('Existing external data bases were found.')
        for ref in xrefs:
            xref_key, xref_val = ref.split(':')
            if not xref_val in id_list:
                if not self.external_sources.has_key(xref_key):
                    self._logger1.info('Creating new external data base node: %s.' % xref_key)
                    new_source, = self.db_connection.data_base.create({'name': xref_key})
                    new_source.add_labels('DB')
                    self.external_sources[xref_key] = new_source
                xref, evidence, link_to = self.db_connection.data_base.create(node({'id': xref_val}),
                                                                     rel(feature_node, 'EVIDENCE', 0),
                                                                     rel(0, 'LINK_TO', self.external_sources[xref_key]))
                xref.add_labels('XRef')
        self._logger1.info('XRefs were created.')

    def search_node(self, label, keys, values):
        if not isinstance(label, str):
            raise ValueError('Label argument must be string.')
        if not keys.__class__.__name__ in 'list':
            raise ValueError('Keys argument must be list.')
        if not values.__class__.__name__ in 'list':
            raise ValueError('Values argument must be list.')
        if not len(keys) == len(values):
            raise TypeError('The lengths of keys and values must be the same.')

        self._logger1.info('Searching for node:')
        session = cypher.Session(self.db_connection.db_link)
        transaction = session.create_transaction()
        query = 'MATCH (node:%s) WHERE ' % label
        for key, value in zip(keys, values):
            if isinstance(value, int) or isinstance(value, float):
                query += 'node.%s=%s AND ' % (key, value)
            else:
                query += 'node.%s="%s" AND ' % (key, value)
        query = query[:-4] + ' RETURN node'
        self._logger1.info(query)
        transaction.append(query)
        try:
            res = list(transaction.commit())[0]
            self._logger1.info('Node was found.')
            return res
        except:
            self._logger1.warning('Node was not found')
            raise UserWarning('Transaction failed.')

class GenomeRelations():

    def __init__(self, db_connection, logger_level = logging.INFO):
        self.db_connection = db_connection
        logging.basicConfig(filename = 'BiomeDB.log', level = logger_level,
                            format = '%(asctime)s %(message)s - %(module)s.%(funcName)s',
                            datefmt='%H:%M:%S-%d.%m.%y')
        self._logger1 = logging.getLogger(__name__)
        self.organism_list = [None, None, None]

    def find_organisms(self):
        return self.db_connection.data_base.find('Organism')

    def set_organism(self, organism):
        if isinstance(organism, str):
            out = self._search_by_name(organism)
        elif isinstance(organism, int):
            out = self._search_by_id(organism)
        else:
            self._logger1.error('Argument must be a string or int.')
            raise ValueError('Argument must be a string or int.'
                             'If string it must contain the name of the organism.'
                             'If int it must be a link to the node in database.')
        self.organism_list = [out, out.get_properties()['name'], node2link(out)]

    def _search_by_name(self, organism):
        session = cypher.Session(self.db_connection.db_link)
        transaction = session.create_transaction()
        query = 'MATCH (org:Organism)' \
                'WHERE org.name="%s" ' \
                'RETURN org' \
                % (organism)
        transaction.append(query)
        transaction_res = list(transaction.commit())[0]
        if not transaction_res:
            err_message = 'There is no organism with name %s in the base.' % organism
            self._logger1.error(err_message)
            raise UserWarning(err_message)
        elif len(transaction_res) > 1:
            err_message = 'There are duplicate organisms with name %s in the base.' % organism
            self._logger1.error(err_message)
            raise UserWarning(err_message)
        else:
            return transaction_res[0][0]

    def _search_by_id(self, organism):
        session = cypher.Session(self.db_connection.db_link)
        transaction = session.create_transaction()
        query = 'START org=node(%d) ' \
                'RETURN org' \
                % (organism)
        transaction.append(query)
        try:
            return list(transaction.commit())[0][0][0]
        except:
            err_message = 'There is no node with link %d in the base.' % organism
            self._logger1.error(err_message)
            raise UserWarning(err_message)

    def _next_overlap_test(self, type_of_node):
        if not isinstance(self.organism_list[0], neo4j.Node):
            self._logger1.error('Organism has not been chosen. First chose the organism by set_organism method.')
            raise ValueError('Organism has not been chosen.')

        if not isinstance(type_of_node, str):
            self._logger1.error('Wrong type of argument.')
            raise ValueError('type_of_node must be string.')

        session = cypher.Session(self.db_connection.db_link)
        transaction = session.create_transaction()
        query = 'START org=node(%s)' \
                'MATCH (org:Organism)<-[:PART_OF]-(ccp) ' \
                'WHERE "Chromosome" IN labels(ccp) OR "Contig" IN labels(ccp) OR "Plasmid" IN labels(ccp) ' \
                'RETURN ccp' \
                % (self.organism_list[2])
        transaction.append(query)
        transaction_res = list(transaction.commit())[0]
        if not transaction_res:
            raise UserWarning('There is no chromosome, contig or plasmid (ccp)'
                              ' for organism "%s" in the base. '
                              'First upload organism to the base '
                              'before using organism.make_next_relation.'
                              % self.organism_list[1])

        check_node_label = self.db_connection.data_base.find(type_of_node)

        try:
            check_node_label = check_node_label.next()
            props = check_node_label.get_properties()
        except:
            self._logger1.error('Label was not found in the base.')
            raise ValueError('Label "%s" is not found in base.' % type_of_node)
        if not props.has_key('start'):
            err_message = 'Label "%s" has not got properties start, end and strand.' % type_of_node
            self._logger1.error(err_message)
            raise ValueError(err_message)

        return transaction_res

    def make_next_relation(self, type_of_node='Feature'):
        """
        The function creates relationships 'NEXT' in the 'data_base' between nodes with label 'type_of_node'.
        The nodes must have property 'start' as the function compare 2 'start' values
        to make the relationship.
        """

        # Make tests and find all ccps for current organism.
        ccps = self._next_overlap_test(type_of_node)
        for ccp in ccps:
            for strand in ('forward', 'reverse'):
                batch = neo4j.WriteBatch(self.db_connection.data_base)
                features = self._feature_start_ordering(ccp[0], strand, type_of_node)
                for i in xrange(1, len(features)):
                    batch.create(rel(features[i-1], 'NEXT', features[i]))
                if ccp[0].get_properties()['type'] == 'circular':
                    try:
                        batch.create(rel(features[-1], 'NEXT', features[0]))
                    except:
                        pass
                batch.submit()
                log_message = 'Created %d relations for %s in %s strand in %s' \
                              % (len(features), type_of_node, strand, self.organism_list[1])
                self._logger1.info(log_message)

    def _feature_start_ordering(self, contig_ref, strand, type_of_node):
        """
        Private method valuecheck suppose to be done by caller method. Never invoke directly!
        """

        contig_ref = node2link(contig_ref)
        session = cypher.Session(self.db_connection.db_link)
        transaction = session.create_transaction()
        if strand == None:
            # Query for OVERLAP
            query = 'START ccp = node(%s) ' \
                    'MATCH (element:%s)-[:PART_OF]->(ccp) ' \
                    'RETURN element ORDER BY element.start' \
                    % (contig_ref, type_of_node)
        else:
            # Query for NEXT
            query = 'START ccp = node(%s) ' \
                    'MATCH (element:%s)-[:PART_OF]->(ccp) ' \
                    'WHERE element.strand="%s" ' \
                    'RETURN element ORDER BY element.start' \
                    % (contig_ref, type_of_node, strand)
        transaction.append(query)
        transaction_out = transaction.commit()
        sorted_features = [result.values[0] for result in transaction_out[0]]
        log_message = 'Amount of features to NEXT/OVERLAP %d' % len(sorted_features)
        self._logger1.info(log_message)
        return sorted_features

    def make_overlap_relation(self, type_of_node='Feature'):
        """
        The function creates relationships 'OVERLAP' in the 'data_base' between nodes with label 'type_of_node'.
        The nodes must have property 'start' as the function compare 2 'start' values
        to make the relationship.
        """

        # Make tests and find all ccps for current organism.
        ccps = self._next_overlap_test(type_of_node)
        for ccp in ccps:
            features = self._feature_start_ordering(ccp[0], None, type_of_node)
            left_edge = [0]
            batch = neo4j.WriteBatch(self.db_connection.data_base)

            # Sorting for OVERLAP
            for i in xrange(2, len(features)):
                for j in xrange(i-1, max(left_edge)-1, -1):
                    if features[j]['end'] >= features[i]['start']:
                        batch.create(rel(features[i], 'OVERLAP', features[j]))
                        # self.db_connection.data_base.create(rel(features[i], 'OVERLAP', features[j]))
                    else:
                        left_edge.append(j)
        batch.submit()
        log_message = 'Created %d relations for %s in %s' \
                      % (len(features), type_of_node, self.organism_list[1])
        self._logger1.info(log_message)

    def connect_organisms_to_taxons(self):
        session = cypher.Session(self.db_connection.db_link)
        transaction = session.create_transaction()
        query = 'MATCH (t:Taxon), (o:Organism) ' \
                'WHERE NOT (t)-[:IS_A]-(o)' \
                'RETURN DISTINCT o.name'
        transaction.append(query)
        transaction_out = transaction.commit()
        if transaction_out:
            for org_name in transaction_out[0]:
                session = cypher.Session(self.db_connection.db_link)
                transaction = session.create_transaction()
                query = 'MATCH (t:Taxon), (o:Organism{name:"%s"}) ' \
                        'WHERE t.scientific_name=o.name ' \
                        'CREATE o-[r:IS_A]->t'\
                        % (org_name[0])
                transaction.append(query)
                try:
                    transaction.commit()
                except:
                    log_message = 'Transaction failed: %s' % query
                    self._logger1.error(log_message)
        else:
            log_message = 'All organisms are connected with taxons.'
            self._logger1.info(log_message)
            print log_message

import doctest
doctest.testfile('test_genbank.txt')