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
    non_gene_list = ['mobile_element', 'repeat_region', 'rep_origin']
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
        self.external_sources = {}

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
                self.create_gene(i)
            elif feature.type in self.non_gene_list or feature.type[:4] == 'misc':
                pass
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
            print 'Pattern found.'
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
                                                                               rel(current_ccp, 'LINK_TO', 0))
                refseq_node.add_labels('XRef')
                print 'XRef was created.'
            else:
                # check relation
                link_to, = self.db_connection.data_base.create(rel(current_ccp, 'LINK_TO', check_xref[0][0]))
                print 'LINK was created.'
        else:
            print 'XRef was found.'


    def feature2node(self):
        if not isinstance(self.rec, SeqRecord.SeqRecord):
            raise ValueError('rec must be Bio.SeqRecord.SeqRecord type.')

        for i in xrange(self.rec.features):
            if self.rec.features[i] == 'gene':
                self.create_gene(self.rec, i)
            method_name = 'self.create' + self.rec.features[i].lower() + '(rec, %d)' % (i)
            try:
                eval(method_name)
            except:
                print self.rec.features[i].type


    def create_gene(self, i):
        # Take info out of gene
        gene = self.rec.features[i]
        gene_dict = {'source': 'GenBank'}
        try:
            gene_name = gene.qualifiers['gene']
            gene_dict['name'] = gene_name
        except:
            gene_name = gene.qualifiers['locus_tag']
            gene_dict['name'] = gene_name

        # If there is no product of gene: RNA or CDS
        if not self.rec.features[i+1].type in ('RNA', 'CDS'):
            try:
                xrefs = gene.qualifiers['db_xref']
            except:
                xrefs = []
            try:
                gene_dict['product'] = gene.qualifiers['product']
            except:
                gene_dict['product'] = []
            # raise UserWarning('gene without product.')

        # Take info out if CDS/RNA
        else:
            if 'CDS' in self.rec.features[i+1].type:
                gene_product = 'cds'
            else:
                gene_product = 'rna'
            gene = self.rec.feaures[i+1]
            gene_dict['locus_tag'] = gene.qualifiers['locus_tag']
            gene_dict['product'] = gene.qualifiers['product']
            try:
                gene_dict['comment'] = gene.qualifiers['note']
            except:
                pass
            try:
                xrefs = gene.qualifiers['db_xref']
            except:
                xrefs = []
            try:
                gene_name = gene.qualifiers['gene']
                gene_dict['name'] = gene_name
            except:
                gene_name = gene.qualifiers['locus_tag']
                gene_dict['name'] = gene_name
        # Find the gene
        gene_dict['start'] = int(gene.location.start) + 1
        gene_dict['end'] = int(gene.location.end)
        gene_dict['strand'] = num2strand(gene.location.strand)
        check_gene = self.search_gene_pattern(gene_dict['start'],
                                              gene_dict['end'],
                                              gene_dict['strand'])
        # If gene is not found
        if not check_gene:
            # Create gene
            gene_node, part_of_org, part_of_ccp =\
                self.db_connection.data_base.create(
                    node(gene_dict),
                    rel(0, 'PART_OF', self.organism_list[0]),
                    rel(0, 'PART_OF', self.ccp_list[0]))
            gene_node.add_labels('Gene', 'Feature', 'DNA', 'BioEntity')
            self.create_xref(xrefs, gene_node)
            try:
                self.create_term(gene_name, gene_node)
            except:
                pass
        else:
            # Update gene
            gene = check_gene[0][0]
            if not gene.get_properties()['name'] == gene_name:
                self.create_term(gene_name, gene)
            source = [gene.get_properties()['source']]
            if not 'GenBank' in source:
                source.append('GenBank')
                gene.update_properties({'source': source})
            self.create_or_update_xref(gene, xrefs)
        invoke = 'self.create_or_update_%s(rec.features[%d], gene)' % (gene_product, i)
        eval(invoke)

    def create_or_update_xref(self, gene, xrefs):
        print 'To be implemented'

    def create_or_update_rna(self, feature, gene):
        print 'To be implemented'

    def create_or_update_cds(self, feature, gene):
        print 'To be implemented'

    def search_gene_pattern(self, start, end, strand):
        if not isinstance(start, int) or not isinstance(end, int):
            raise ValueError('start and end must be string.')
        if not isinstance(strand, str):
            raise ValueError('strand must be string.')
        if not start in ["linear", "circular", "unknown"]:
            raise ValueError('strand must be "linear", "circular" or "unknown".')

        session = cypher.Session(self.db_connection.db_link)
        transaction = session.create_transaction()
        query = 'START org=node(%s), ccp=node(%s) ' \
                'MATCH (org)<-[:PART_OF]-(ccp)<-[:PART_OF]-(g:Gene)-[:PART_OF]->(org) ' \
                'WHERE g.start=%d, g.end=%d, g.strand=%d ' \
                'RETURN g' \
                % (self.organism_list[2], self.ccp_list[2], start, end, strand)
        transaction.append(query)
        try:
            transaction_res = list(transaction.commit())[0]
            return transaction_res
        except:
            print 'Transaction failed.'

    def create_term(self, name, bioentity):
        if not isinstance(name, str):
            raise ValueError('id must be string.')
        if not isinstance(bioentity, neo4j.Node):
            raise ValueError('bioentity must be py2neo.neo4j.Node.')

        term, has_name = self.db_connection.data_base.create(node({'text': name}),
                                                             rel(bioentity, 'HAS_NAME', 0))
        term.add_labels('Term')
        # self.db_connection.data_base.create(rel(bioentity, 'EVIDENCE', xref))

    def create_xref(self, refs, feature):
        """
        Method creates XRef to a feature with relation EVIDENCE.
        :param refs: list
        :param feature: py2neo.neo4j.Node
        :return: None
        """
        if not isinstance(refs, list):
            raise ValueError('refs must be a list.')
        if not isinstance(feature, neo4j.Node):
            raise ValueError('feature must be py2neo.neo4j.Node.')

        for ref in refs:
            xref_key, xref_val = ref.split(':')
            if not self.external_sources.has_key(xref_key):
                new_source, = self.db_connection.data_base.create({'name': xref_key})
                new_source.add_label('DB')
                self.external_sources[xref_key] = new_source
            xref, evidence = self.db_connection.data_base.create(node({'id': xref_val}),
                                                                 rel(feature, 'EVIDENCE', 0),
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

import doctest
doctest.testfile('test_genbank.txt')