#from ...api import *
from py2neo import neo4j, node, rel, cypher
import biome.load.genbank.genbank as gb
import os
import warnings

def update_source_property(node):
    if not isinstance(node, gb.neo4j.Node):
        raise TypeError('The node argument must be an object of neo4j.Node class!')
    source = node.get_properties()['source']
    if 'RegulonDB' in source:
        pass
    elif isinstance(source, basestring):
        node.update_properties({'source': [source, 'RegulonDB']})
    elif isinstance(source, list):
        node.update_properties({'source': source.append('RegulonDB')})
    else:
        raise Exception('Unexpected source type!')


class RegulonDB():
    """

    """
    def __init__(self, directory, ecoli_name='Escherichia coli str. K-12 substr. MG1655',
                 chro_name='Escherichia coli str. K-12 substr. MG1655, complete genome.',
                 dblink='http://localhost:7474/db/data/'):
        if not isinstance(ecoli_name, basestring):
            raise TypeError('The ecoli_name argument must be a string!')
        if not isinstance(dblink, basestring):
            raise TypeError('The connection argument must be a string!')
        if not os.path.isdir(directory):
            raise ValueError('The directory does not exist!')
        self.directory = directory
        self.ecoli_name = ecoli_name
        self.chro_name = chro_name
        self.dblink = dblink
        self.connection = neo4j.GraphDatabaseService(self.dblink)

        try:
            ecoli_node = list(self.connection.find('Organism', 'name', self.ecoli_name))
        except:
            raise ValueError('Check the dblink! Could not connect!')

        if not ecoli_node:
            raise ValueError('There is no organism node with %s name!' % self.ecoli_name)
        self.ecoli_node = ecoli_node[0]

        try:
            chro_node = list(self.connection.find('Chromosome', 'name', self.chro_name))
        except:
            raise ValueError('Check the dblink! Could not connect!')

        if not chro_node:
            raise ValueError('There is no chromosome node with %s name!' % self.chro_name)
        self.chro_node = chro_node[0]

    def __repr__(self):
        return "RegulonDB object for %s\nLink to database: %s" \
               % (self.ecoli_name, self.dblink)

    def __str__(self):
        return "RegulonDB object for %s\nLink to database: %s" \
               % (self.ecoli_name, self.dblink)

    def check_create_terms(self, bioentity, name):
        if not isinstance(bioentity, gb.neo4j.Node):
            raise TypeError('The node argument must be an object of neo4j.Node class!')
        if bioentity['name'] != name:
            term, rel_pro = self.connection.create(
                node({'text': name}),
                rel(0, 'HAS_NAME', bioentity))
            term.add_labels('Term')

    def create_operons(self):
        f = open(self.directory + 'Operons.txt', 'r')
        data = f.readlines()
        f.close()
        i = 0
        for line in data:
            if line[0] == '#':
                continue
            chunks = line.split('\t')
            operon, term, term_rel, org_rel = self.connection.\
                create(node({'name': chunks[0], 'start': chunks[1],
                             'end': chunks[2], 'strand': chunks[3],
                             'evidence': chunks[6], 'source': 'RegulonDB'}),
                       node({'text': chunks[0]}),
                       rel(0, 'HAS_NAME', 1),
                       rel(0, 'PART_OF', self.ecoli_node))
            operon.add_labels('Operon', 'BioEntity', 'DNA')
            i += 1
        print '%d operons were created!' % i

    def create_update_promoters(self):
        f = open(self.directory + 'All Promoters.txt', 'r')
        data = f.readlines()
        f.close()
        notfound = 0
        updated = 0
        for line in data:
            if line[0] == '#':
                continue
            regid, name, strand, tss, sigma, seq, evidence = line.split('\t')
            tss = int(tss)

            # skipping incomplete data
            if '' in [regid, name, strand, tss]:
                continue

            query = 'MATCH (ch:Chromosome {name: "%s"})<-[:PART_OF]-' \
                    '(p:Promoter {tss: %d})-[:PART_OF]->' \
                    '(o:Organism {name: "%s"}) ' \
                    'RETURN p' % (self.chro_name, tss,  self.ecoli_name)
            promoter_query = neo4j.CypherQuery(self.connection, query)
            promoter_nodes = promoter_query.execute()

            # no promoter with the tss
            if not promoter_nodes:
                promoter, term, rel_org, rel_chr, rel_term = self.connection.create(
                    node({'name': name, 'start': tss,
                          'end': tss, 'strand': strand,
                          'tss': tss, 'seq': seq,
                          'evidence': evidence, 'Reg_id': regid,
                          'source': 'RegulonDB'}),
                    node({'text': name}),
                    rel(0, 'PART_OF', self.ecoli_node),
                    rel(0, 'PART_OF', self.chro_node),
                    rel(0, 'HAS_NAME', 1))
                promoter.add_labels('Promoter', 'Feature', 'BioEntity', 'DNA')
                term.add_labels('Term')
                notfound += 1
            else:
                # one promoter with the tss
                for record in promoter_nodes.data:
                    promoter = record.values[0]
                    promoter.update_properties({'seq': seq,
                                                'evidence': evidence,
                                                'Reg_id': regid})
                    update_source_property(promoter)
                    self.check_create_terms(promoter, name)
                    updated += 1

                # duplicates!
                if len(promoter_nodes.data) > 1:
                    warnings.warn("There are %d nodes for a promoter with tss"
                                  " in the %d position! It was skipped!\n"
                                  % (len(promoter_nodes.data), tss))

        print '%d promoters were updated!\n' \
              '%d promoters were created!' % (updated, notfound)




