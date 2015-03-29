#from ...api import *
from py2neo import neo4j, node, rel, cypher
import biome.load.genbank.genbank as gb
import os

class RegulonDB():
    """

    """
    def __init__(self, directory, ecoli_name='Escherichia coli str. K-12 substr. MG1655',
                 dblink='http://localhost:7474/db/data/'):
        if not isinstance(ecoli_name, basestring):
            raise TypeError('The ecoli_name argument must be a string!')
        if not isinstance(dblink, basestring):
            raise TypeError('The connection argument must be a string!')
        if not os.path.isdir(directory):
            raise ValueError('The directory does not exist!')
        self.directory = directory
        self.ecoli_name = ecoli_name
        self.dblink = dblink
        self.connection = neo4j.GraphDatabaseService(self.dblink)

        try:
            ecoli_node = list(self.connection.find('Organism', 'name', self.ecoli_name))
        except:
            raise ValueError('Check the dblink! Could not connect!')

        if not ecoli_node:
            raise ValueError('There is no organism node with %s name!' % self.ecoli_name)

        self.ecoli_node = ecoli_node[0]

    def __repr__(self):
        return "RegulonDB object for %s\nLink to database: %s" \
               % (self.ecoli_name, self.dblink)

    def __str__(self):
        return "RegulonDB object for %s\nLink to database: %s" \
               % (self.ecoli_name, self.dblink)

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
                             'evidence': chunks[6]}),
                       node({'text': chunks[0]}),
                       rel(0, 'HAS_NAME', 1),
                       rel(0, 'PART_OF', self.ecoli_node))
            operon.add_labels('Operon', 'BioEntity', 'DNA')
            i += 1
        print '%d operons were created!' % i

    def create_update_promoters(self):
        f = open(self.directory + 'Transcription Units.txt', 'r')
        data = f.readlines()
        f.close()
        notfound = 0
        i = 0
        for line in data:
            if line[0] == '#':
                continue
            chunks = line.split('\t')

    def create_update_tus(self):
        f = open(self.directory + 'Transcription Units.txt', 'r')
        data = f.readlines()
        f.close()
        notfound = 0
        i = 0
        for line in data:
            if line[0] == '#':
                continue
            chunks = line.split('\t')
            tu_nodes = list(self.connection.find('TU', 'name', chunks[1]))

            if not tu_nodes:
                tu = self.connection.create(
                    node({'name': chunks[0], 'start': chunks[1],
                          'end': chunks[2], 'strand': chunks[3],
                          'evidence': chunks[6]}))
            elif len(tu_nodes)==1:
                tu = tu_nodes[0]
            else:
                print "%d TUs with the same name were found" % len(tu_nodes)







