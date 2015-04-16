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


def tf_effect(effect):
    if effect == '+':
        return 'ACTIVATES'
    if effect == '-':
        return 'REPRESSES'
    if effect == '+-':
        return 'DUAL'
    if effect == '?':
        return 'UNKNOWN'


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

    def relation_with_tu(self, tu_name, element):
        query = 'MATCH (o:Organism {name: "%s"})<-[:PART_OF]-' \
                    '(tu:TU)-[:HAS_NAME]->(:Term {text: "%s"}) ' \
                    'RETURN tu' % (self.ecoli_name, tu_name)
        res = neo4j.CypherQuery(self.connection, query)
        res_nodes = res.execute()

        if not res_nodes:
            warnings.warn("There is no node for a TU with name %s!\n"
                              "It was skipped!\n" % tu_name)
            return 1
        else:
            for tu in res_nodes.data:
                rel_tu = self.connection.create(
                    rel(tu.values[0], 'CONTAINS', element))
            return 0

    def create_operons(self):
        f = open(self.directory + 'Operons.txt', 'r')
        data = f.readlines()
        f.close()
        i = 0
        for line in data:
            if line[0] == '#':
                continue
            chunks = line.split('\t')

            ### testing
            if chunks[0] == '' or chunks[1] == '' or chunks[2] == 0:
                continue
            if chunks[3] == '':
                chunks[3] = 'unknown'

            operon, term, term_rel, org_rel = self.connection.\
                create(node({'name': chunks[0], 'start': int(chunks[1]),
                             'end': int(chunks[2]), 'strand': chunks[3],
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
        created, updated = [0]*2

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
            res = neo4j.CypherQuery(self.connection, query)
            res_nodes = res.execute()

            # creating promoter
            if not res_nodes:
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
                created += 1
            else:
                # one promoter with the tss
                for record in res_nodes.data:
                    promoter = record.values[0]
                    promoter.update_properties({'seq': seq,
                                                'evidence': evidence,
                                                'Reg_id': regid})
                    update_source_property(promoter)
                    self.check_create_terms(promoter, name)
                    updated += 1

                # duplicates!
                if len(res_nodes.data) > 1:
                    warnings.warn("There are %d nodes for a promoter with tss"
                                  " in the %d position! It was skipped!\n"
                                  % (len(res_nodes.data), tss))

        print '%d promoters were updated!\n' \
              '%d promoters were created!' % (updated, created)

    def create_update_tus(self):
        f = open(self.directory + 'Transcription Units.txt', 'r')
        data = f.readlines()
        f.close()
        created, updated, problem = [0]*3
        for line in data:
            if line[0] == '#':
                continue
            regid, name, operon, genes_name, pro, evidence = line.split('\t')

            ### testing
            if '' in [regid, operon]:
                continue

            # searching for TU with the same name
            query = 'MATCH (t:Term {text: "%s"})<-[:HAS_NAME]-' \
                    '(p:Promoter)<-[:CONTAINS]-(tu:TU)-[:PART_OF]->' \
                    '(o:Organism {name: "%s"}) ' \
                    'RETURN tu' % (pro, self.ecoli_name)
            res = neo4j.CypherQuery(self.connection, query)
            res_nodes = res.execute()

            # no tu with the name was found
            if not res_nodes:
                tu, term, rel_org, rel_term = self.connection.create(
                    node({'name': name, 'evidence': evidence,
                          'Reg_id': regid, 'source': 'RegulonDB'}),
                    node({'text': name}),
                    rel(0, 'PART_OF', self.ecoli_node),
                    rel(0, 'HAS_NAME', 1))
                tu.add_labels('TU', 'BioEntity', 'DNA')
                term.add_labels('Term')
                created += 1

                # creating a relation (:TU)-[:CONTAINS]->(:Promoter)
                query = 'MATCH (t:Term {text: "%s"})<-[:HAS_NAME]-' \
                        '(p:Promoter)-[:PART_OF]->(o:Organism {name: "%s"}) ' \
                        'RETURN p' % (pro, self.ecoli_name)
                res = neo4j.CypherQuery(self.connection, query)
                res_nodes = res.execute()

                if not res_nodes:
                    warnings.warn("There is no node for a promoter with name %s!\n"
                                  "It was skipped!\n" % pro)

                # if there are promoters-duplicates
                elif len(res_nodes) > 1:
                    warnings.warn("There are %d nodes for a promoter with name %s!\n"
                                  "They were skipped!\n" % (len(res_nodes), pro))
                else:
                    rel_promoter = self.connection.create(
                        rel(tu, 'CONTAINS', res_nodes.data[0].values[0]))

            elif len(res_nodes.data) == 1:
                tu = res_nodes.data[0].values[0]
                tu.update_properties({'evidence': evidence,
                                      'Reg_id': regid})
                update_source_property(tu)
                self.check_create_terms(tu, name)
                updated += 1
            else:
                problem += 1
                warnings.warn("There are %d nodes for a TU with name %s! All of them have the same promoter %s!\n"
                              "They were skipped!\n" % (len(res_nodes.data), name, pro))
                continue

            # creating a relation (:TU)<-[:CONTAINS]-(:Operon)
            operon_node = list(self.connection.find('Operon', 'name', operon))

            if not operon_node:
                warnings.warn("There is no node for an operon with name %s!\n"
                              "It was skipped!\n" % operon)

            # if there are operons-duplicates
            elif len(operon_node) > 1:
                warnings.warn("There are %d nodes for an operon with name %s!\n"
                              "They were skipped!\n" % (len(operon_node), operon))
            else:
                rel_operon = self.connection.create(
                    rel(operon_node[0], 'CONTAINS', tu))


        print "%d TUs were updated and connected to operons!\n" \
              "%d TUs were created and connected to operons!\n" \
              "There were problems with %d TUs." % (updated, created, problem)

    def create_update_terminators(self):
        f = open(self.directory + 'Terminators.txt', 'r')
        data = f.readlines()
        f.close()
        created, updated, problem = [0]*3
        for line in data:
            if line[0] == '#':
                continue
            regid, start, end, strand, seq, tu, type, operon, ref, evidence = line.split('\t')
            start, end = [int(start), int(end)]

            # skipping incomplete data
            if '' in [regid, strand, start, end] or 0 in [start, end]:
                continue

            query = 'MATCH (ch:Chromosome {name: "%s"})<-[:PART_OF]-' \
                    '(t:Terminator {start: %d, end: %d, strand: "%s"}) ' \
                    'RETURN t' % (self.chro_name, start, end, strand)
            res = neo4j.CypherQuery(self.connection, query)
            res_nodes = res.execute()

            # creating terminator
            if not res_nodes:
                terminator, rel_chr = self.connection.create(
                    node({'start': start, 'end': end,
                          'strand': strand, 'seq': seq,
                          'evidence': evidence, 'Reg_id': regid,
                          'source': 'RegulonDB'}),
                    rel(0, 'PART_OF', self.chro_node))
                terminator.add_labels('Terminator', 'Feature', 'DNA')
                created += 1

            elif len(res_nodes.data) == 1:
                    terminator = res_nodes.data[0].values[0]
                    terminator.update_properties({'seq': seq,
                                                  'evidence': evidence,
                                                  'Reg_id': regid})
                    update_source_property(terminator)
                    updated += 1

            # duplicates!
            else:
                warnings.warn("There are %d nodes for a terminator with "
                              "location (%d, %d, %s)! It was skipped!\n"
                              % (len(res_nodes.data), start, end, strand))
                continue

            # creating relations (:TU)-[:CONTAINS]->(:Terminator)
            rel_tu = self.relation_with_tu(tu, terminator)
            problem = problem + rel_tu

        print '%d terminators were updated!\n' \
              '%d terminators were created!\n' \
              'There were problems with %d terminators.' \
              % (updated, created, problem)

    def create_update_BSs(self):
        f = open(self.directory + 'TF binding sites.txt', 'r')
        data = f.readlines()
        f.close()
        created, updated, problem = [0]*3

        for line in data:
            if line[0] == '#':
                continue
            regid, name, site_id, start, end, strand, inter_id, tu_name, effect, pro, center, seq, evidence = line.split('\t')

            ### testing
            if '' in [regid, strand, start, end, center] or 0 in [start, end]:
                continue

            start, end, center = [int(start), int(end), float(center)]

            query = 'MATCH (o:Organism {name: "%s"})<-[:PART_OF]-' \
                    '(tu:TU)-[:HAS_NAME]-(t1:Term {text: "%s"}), ' \
                    '(tu)-[:CONTAINS]->(p:Promoter)-[:HAS_NAME]-(t2:Term {text: "%s"}) ' \
                    'RETURN p, tu' % (self.ecoli_name, tu_name, pro)
            res = neo4j.CypherQuery(self.connection, query)
            res_nodes = res.execute()

            if not res_nodes:
                problem += 1
                continue
            elif len(res_nodes.data) == 1:
                promoter = res_nodes.data[0].values[0]
                tu = res_nodes.data[0].values[1]
            else:
                warnings.warn("It is impossible to identify a transcription "
                              "unit for a binding site with location "
                              "(%d, %d, %s)! It was skipped!\n"
                              % (start, end, strand))
                continue

            ### calculating BS position in MetaCyc
            MC_start = int(promoter["tss"] + center)

            query = 'MATCH (o:Organism {name: "%s"})<-[:PART_OF]-' \
                    '(tu:TU)-[:HAS_NAME]-(t1:Term {text: "%s"}), ' \
                    '(tu)-[:CONTAINS]->(p:Promoter)-[:HAS_NAME]-(t2:Term {text: "%s"}), ' \
                    '(tu)-[:CONTAINS]->(bs:BS {start: %d, strand: "%s"}) ' \
                    'RETURN bs' % (self.ecoli_name, tu_name, pro, MC_start, strand)
            res = neo4j.CypherQuery(self.connection, query)
            res_nodes = res.execute()

            print query
            print res_nodes.data

            # creating BS
            if not res_nodes:
                bs, transreg, rel_bs_transreg, rel_chr, rel_pro, rel_tu = self.connection.create(
                    node({'start': start, 'end': end,
                          'strand': strand, 'seq': seq,
                          'evidence': evidence, 'Reg_id': site_id,
                          'source': 'RegulonDB', 'center': center}),
                    node({'Reg_id': inter_id, 'source': 'RegulonDB'}),
                    rel(0, 'PARTICIPATES_IN', 1),
                    rel(0, 'PART_OF', self.chro_node),
                    rel(1, tf_effect(effect), promoter),
                    rel(tu, 'CONTAINS', 0))
                bs.add_labels('BS', 'Feature', 'DNA')
                transreg.add_labels('TranscriptionRegulation', 'RegulationEvent', 'Binding')
                created += 1

            else:
                for record in res_nodes.data:
                    bs = record.values[0]
                    bs.update_properties({'seq': seq, 'start': start,
                                          'end': end, 'strand': strand,
                                          'evidence': evidence,
                                          'Reg_id': site_id,
                                          'center': center})
                    update_source_property(bs)

                    transreg, rel_bs_transreg, rel_pro = self.connection.create(
                        node({'Reg_id': inter_id, 'source': 'RegulonDB'}),
                        rel(bs, 'PARTICIPATES_IN', 0),
                        rel(0, tf_effect(effect), promoter))
                    transreg.add_labels('TranscriptionRegulation', 'RegulationEvent', 'Binding')
                    updated += 1

            # creating relations (:TF)-[:PARTICIPATES_IN]->(:TranscriptionRegulation)

        print '%d BSs were updated!\n' \
              '%d BS were created!\n' \
              'There were problems with %d BSs.' \
              % (updated, created, problem)