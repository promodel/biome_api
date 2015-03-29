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
            # session = cypher.Session(self.dblink)
            # transaction = session.create_transaction()
            # query = 'MATCH (org:Organism {name:"%s"}) ' \
            #         'RETURN org' % self.ecoli_name
            # transaction.append(query)
            # transaction_res = transaction.commit()[0]
            ecoli_node = list(self.connection.find('Organism', 'name', self.ecoli_name))
            # print transaction_res[0][0]
        except:
            raise ValueError('Check the dblink! Could not connect!')

        if not ecoli_node:
            raise ValueError('There is no organism node with %s name!' % self.ecoli_name)

        self.ecoli_node = ecoli_node


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
        for line in data:
            if line[0] == '#':
                continue
            chunks = line.split('\t')






