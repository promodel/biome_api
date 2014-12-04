from Bio import SeqIO, SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from py2neo import neo4j, node, rel, cypher
from time import ctime, time
# import datetime
import warnings
import logging
import os
# from fileinput import input
# from Bio.GenBank import RecordParser as rp

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
    def __init__(self, db_link='http://localhost:7474/db/data/'):
        if not isinstance(db_link, str):
            raise TypeError ('Argument must be string.')

        self.db_link = db_link
        self.data_base = neo4j.GraphDatabaseService(self.db_link)

class GenBank():
    non_gene_list = ['mobile_element', 'repeat_region', 'rep_origin', 'STS']
    gene_product_list = ['CDS', 'tRNA', 'rRNA', 'tmRNA', 'ncRNA', 'mRNA']
    def __init__(self, gb_file, db_connection):
        if not isinstance(gb_file, str):
            raise TypeError('gb_file must be a string with filename.')
        if not isinstance(db_connection, BioGraphConnection):
            raise TypeError('db_connection must be an instance of the class BioGraphConnection.')
        if not os.path.isfile(gb_file):
            raise ValueError('There is no %s in current directory.' % gb_file)

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
        for i, feature in enumerate(self.rec.features):
            if feature.type == 'gene':
                self.make_gene_and_product(i)
            elif feature.type in self.non_gene_list or feature.type[:4] == 'misc':
                self.make_non_gene(feature)
            elif feature.type in self.gene_product_list:
                pass
            else:
                print 'Unknown element %s was skipped.' % feature.type

    def _read_gb_file(self):
        rec = SeqIO.read(self.gb_file, 'genbank')
        if not rec:
            raise ValueError('gb-file is empty.')
        seq_file = open(self.gb_file, 'r')
        first_line = seq_file.readline()
        seq_file.close()
        if 'circular' in first_line:
            seq_type = 'circular'
        else:
            seq_type = 'linear'
        return rec, seq_type

    def search_pattern_cypher(self):
        taxon = self.rec.features[0].qualifiers['db_xref'][0].split(':')[1]
        session = cypher.Session(self.db_connection.db_link)
        transaction = session.create_transaction()
        query = 'MATCH (org:Organism)<-[:PART_OF]-(ccp:Chromosome)-[:EVIDENCE]->(xref:XRef) ' \
                'WHERE org.name="%s" AND ccp.length=%d AND xref.id="%s" RETURN org, ccp, xref' \
                % (self.rec.annotations['source'],
                   int(self.rec.features[0].location.end),
                   self.rec.annotations['accessions'][0])
        transaction.append(query)
        transaction_res = list(transaction.commit())[0]
        if transaction_res:
            transaction_res = transaction_res[0]
            self.organism_list = [transaction_res[0], self.rec.annotations['organism'], node2link(transaction_res[0])]
            self.ccp_list = [transaction_res[1], self.rec.description, node2link(transaction_res[1])]
        else:
            # Create or find organism
            self.create_or_update_organism()

            # Create or find nodes for chromosome, contig or plasmid
            self.create_or_update_ccp()

    def create_or_update_organism(self):
        # searching organism
        organism_name = self.rec.annotations['organism']
        search_organism = list(self.db_connection.data_base.find('Organism', 'name', organism_name))
        if not search_organism:
            # creating organism
            current_organism, = self.db_connection.data_base.create(
                node({'name': organism_name}))
            current_organism.add_labels('Organism')
        elif len(search_organism) > 1:
            raise Exception('There are duplicates in the DB.')
        else:
            current_organism = search_organism[0]
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
            raise UserWarning('Unknown genome element')
        ccp_length = len(self.rec)
        ccp_name = description
        search_ccp = self.search_node(ccp_label, ['name', 'length'],
                                      [ccp_name, ccp_length])
        if not search_ccp:
            # Creating chromosome, contig or plasmid
            current_ccp, = self.db_connection.data_base.create(node({'name': ccp_name,
                                                                   'length': ccp_length,
                                                                   'type': self.seq_type}))
            # Adding label
            current_ccp.add_labels(ccp_label, 'BioEntity', 'DNA')
            print 'Chromosome was created.'
        else:
            current_ccp = search_ccp[0][0]
            print 'Chromosome was found.'
        self.ccp_list = [current_ccp, ccp_name, node2link(current_ccp)]
        self.db_connection.data_base.create(rel(self.ccp_list[0], 'PART_OF', self.organism_list[0]))

        # Check XRefs
        refseq = self.rec.annotations['accessions'][0]
        session = cypher.Session(self.db_connection.db_link)
        transaction = session.create_transaction()
        query = 'START ccp=node(%s) MATCH (ccp)-[e:EVIDENCE]->(xref) WHERE xref.id="%s" RETURN e' \
                % (self.ccp_list[2], refseq)
        transaction.append(query)
        transaction_res = list(transaction.commit())[0]
        if not transaction_res:
            check_xref = self.search_node('XRef', ['id'], [refseq])
            if not check_xref:
                # Create XRefs
                refseq_node, refseq_link = self.db_connection.data_base.create(node({'id': refseq}),
                                                                               rel(current_ccp, 'EVIDENCE', 0))
                refseq_node.add_labels('XRef')
                print 'XRef was created.'
            else:
                # check relation
                link_to, = self.db_connection.data_base.create(rel(current_ccp, 'LINK_TO', check_xref[0][0]))
                print 'LINK was created.'
        else:
            print 'XRef was found.'


    def make_gene_and_product(self, i):
        gene = self.rec.features[i]
        # Take info out of gene
        gene_dict = {'locus_tag':gene.qualifiers['locus_tag'][0]}
        try:
            next_feature = self.rec.features[i+1]
            if not next_feature.type in self.gene_product_list:
                next_feature = self.rec.features[i]
            else:
                gene_dict['product'] = next_feature.qualifiers['product'][0]
                pass
        except:
            next_feature = self.rec.features[i]
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
            self.create_xref(xrefs=xrefs, feature_node=gene_node, check=check)
        except:
            pass

        if next_feature.type == 'CDS':
            self.create_or_update_cds(next_feature, gene_node)
        elif next_feature.type in self.gene_product_list[1:-1]:
            self.create_or_update_rna(next_feature, gene_node)
        else:
            pass


    def create_or_update_gene(self, gene_dict):
        # Searching the gene in the DB
        check_gene = self.search_gene_pattern(gene_dict['start'],
                                              gene_dict['end'],
                                              gene_dict['strand'])
        # If gene is not found
        if not check_gene:
            print 'Creating gene'
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
            print 'Gene is found. Updating.'

            # Updating the gene
            gene_node = check_gene[0][0]
            gene_props = gene_node.get_properties()

            # Updating source
            self.update_source(gene_node, gene_props)

            # Checking locus_tag
            try:
                gene_locus_tag = gene_props['locus_tag']
                if gene_locus_tag != gene_dict['locus_tag']:
                    print 'Locus_tag is different' # Log
            except:
                gene_node.update_properties({'locus_tag': gene_dict['locus_tag']})

            # Checking product
            try:
                gene_product = gene_props['product']
                if gene_product != gene_dict['product']:
                    print 'Product is different' # Log
            except:
                try:
                    gene_node.update_properties({'product': gene_dict['product']})
                except:
                    pass

            # Flag for create_xref method.
            updated = True

            # Checking name of existing gene
            if gene_props['name'] != gene_dict['name']:
                self.create_term(gene_dict['name'], gene_node)
        return gene_node, updated

    def create_or_update_rna(self, feature, gene_node):
        # Taking info for the RNA
        try:
            rna_dict = {'name':feature.qualifiers['product'][0]}
        except:
            rna_dict = {'name':feature.qualifiers['locus_tag'][0]}

        try:
            rna_dict['comment'] = feature.qualifiers['note'][0]
        except:
            pass

        # Searching for the RNA
        session = cypher.Session(self.db_connection.db_link)
        transaction = session.create_transaction()
        query = 'START g=node(%s) ' \
                'MATCH (g)-[:ENCODES]->(r:RNA) ' \
                'RETURN r' % node2link(gene_node)
        transaction.append(query)
        check_rna = transaction.commit()[0]
        if not check_rna:
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
            xrefs = [db for db in feature.qualifiers['db_xref'] if db.split(':')[0] not in ['GeneID']]
            self.create_xref(xrefs=xrefs, feature_node=rna_node, check=updated),
        except:
            pass

    def create_or_update_cds(self, feature, gene_node):
        # Taking info for the polypetide
        poly_dict = {'name':feature.qualifiers['product'][0]}
        try:
            seq = feature.qualifiers['translation'][0]
        except:
            dna_seq = self.rec.seq[feature.location.nofuzzy_start:feature.location.nofuzzy_end]
            seq = dna_seq.translate(feature.qualifiers['trans_table'][0])[:-1]
        poly_dict['seq'] = seq
        try:
            poly_dict['comment'] = feature.qualifiers['note'][0]
        except:
            pass

        # Searching for the polypeptide
        session = cypher.Session(self.db_connection.db_link)
        transaction = session.create_transaction()
        query = 'START g=node(%s) ' \
                'MATCH (g)-[:ENCODES]->(p:Polypeptide) ' \
                'WHERE p.seq="%s" ' \
                'RETURN p' % (node2link(gene_node), seq)
        transaction.append(query)
        check_poly = transaction.commit()[0]
        if not check_poly:
            poly_dict['source'] = ['GenBank']
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
            xrefs = [db for db in feature.qualifiers['db_xref'] if db.split(':')[0] not in ['GeneID']]
            self.create_xref(xrefs=xrefs, feature_node=poly_node, check=updated),
        except:
            pass

    def make_non_gene(self, feature):
        feature_dict ={}
        feature_dict['start'] = int(feature.location.start) + 1
        feature_dict['end'] = int(feature.location.end)
        feature_dict['strand'] = num2strand(feature.location.strand)
        check_feature = self.search_non_gene_pattern(feature_dict['start'],
                                              feature_dict['end'],
                                              feature_dict['strand'])
        try:
            feature_dict['comment'] = feature.qualifiers['note'][0]
        except:
            pass
        try:
            feature_dict['type'] = feature.qualifiers['mobile_element_type']
        except:
            pass
        if not check_feature:
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
            xrefs = feature.qualifiers['db_xref']
            self.create_xref(xrefs=xrefs, feature_node=feature_node, check=updated),
        except:
            pass

    def update_source(self, node, props):
        source = props['source']
        if not isinstance(source, list):
            source = [source]
        if not 'GenBank' in source:
            source.append('GenBank')
            node.update_properties({'source': sorted(source)})

    def search_gene_pattern(self, start, end, strand):
        if not isinstance(start, int) or not isinstance(end, int):
            raise ValueError('start and end must be string.')
        if not isinstance(strand, str):
            raise ValueError('strand must be string.')
        if not strand in ["forward", "reverse", "unknown"]:
            raise ValueError('strand must be "forward", "reverse" or "unknown".')

        session = cypher.Session(self.db_connection.db_link)
        transaction = session.create_transaction()
        query = 'START org=node(%s), ccp=node(%s) ' \
                'MATCH (ccp)<-[:PART_OF]-(g:Gene)-[:PART_OF]->(org) ' \
                'WHERE g.start=%d AND g.end=%d AND g.strand="%s" ' \
                'RETURN g' \
                % (self.organism_list[2], self.ccp_list[2], start, end, strand)
        transaction.append(query)
        try:
            transaction_res = list(transaction.commit())[0]
            return transaction_res
        except:
            print 'Transaction failed.'

    def search_non_gene_pattern(self, start, end, strand):
        if not isinstance(start, int) or not isinstance(end, int):
            raise ValueError('start and end must be string.')
        if not isinstance(strand, str):
            raise ValueError('strand must be string.')
        if not strand in ["forward", "reverse", "unknown"]:
            raise ValueError('strand must be "forward", "reverse" or "unknown".')

        session = cypher.Session(self.db_connection.db_link)
        transaction = session.create_transaction()
        query = 'START ccp=node(%s) ' \
                'MATCH (ccp)<-[:PART_OF]-(ng) ' \
                'WHERE ng.start=%d AND ng.end=%d AND ng.strand="%s" ' \
                'RETURN ng' \
                % (self.organism_list[2], start, end, strand)
        transaction.append(query)
        try:
            transaction_res = list(transaction.commit())[0]
            return transaction_res
        except:
            print 'Transaction failed.'

    def create_term(self, name, bioentity):
        if not isinstance(name, str):
            raise ValueError('name must be string.')
        if not isinstance(bioentity, neo4j.Node):
            raise ValueError('bioentity must be py2neo.neo4j.Node.')

        term, has_name = self.db_connection.data_base.create(node({'text': name}),
                                                             rel(bioentity, 'HAS_NAME', 0))
        term.add_labels('Term')
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

        id_list = []
        if check == True:
            session = cypher.Session(self.db_connection.db_link)
            transaction = session.create_transaction()
            query = 'START feature=node(%s)' \
                    ' MATCH (feature)-[:EVIDENCE]->(xref:XRef)-[:LINK_TO]->(db)' \
                    ' RETURN xref.id, db.name, db' % node2link(feature_node)
            transaction.append(query)
            transaction_res = list(transaction.commit())[0]
            for rec in transaction_res:
                id_list.append(rec[0])
                if not self.external_sources.has_key(rec[1]):
                    self.external_sources[rec[1]] = rec[2]
        for ref in xrefs:
            xref_key, xref_val = ref.split(':')
            if not xref_val in id_list:
                if not self.external_sources.has_key(xref_key):
                    new_source, = self.db_connection.data_base.create({'name': xref_key})
                    new_source.add_labels('DB')
                    self.external_sources[xref_key] = new_source
                xref, evidence, link_to = self.db_connection.data_base.create(node({'id': xref_val}),
                                                                     rel(feature_node, 'EVIDENCE', 0),
                                                                     rel(0, 'LINK_TO', self.external_sources[xref_key]))
                xref.add_labels('XRef')

    def search_node(self, label, keys, values):
        if not isinstance(label, str):
            raise ValueError('Label argument must be string.')
        if not keys.__class__.__name__ in 'list':
            raise ValueError('Keys argument must be list.')
        if not values.__class__.__name__ in 'list':
            raise ValueError('Values argument must be list.')
        if not len(keys) == len(values):
            raise TypeError('The lengths of keys and values must be the same.')

        session = cypher.Session(self.db_connection.db_link)
        transaction = session.create_transaction()
        query = 'MATCH (node:%s) WHERE ' % label
        for key, value in zip(keys, values):
            if isinstance(value, int) or isinstance(value, float):
                query += 'node.%s=%s AND ' % (key, value)
            else:
                query += 'node.%s="%s" AND ' % (key, value)
        query = query[:-4] + ' RETURN node'
        transaction.append(query)
        try:
            return list(transaction.commit())[0]
        except:
            raise UserWarning('Transaction failed.')

# import doctest
# doctest.testfile('test_genbank.txt')