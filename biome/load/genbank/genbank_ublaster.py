import logging
from time import time
from py2neo import neo4j, node, rel, cypher
import warnings
from genbank import BioGraphConnection, node2link
from genbank_blaster import MakeJob
import hashlib

class BlastUploader():
    def __init__(self, db_connection, logger_level=logging.INFO):
        if not isinstance(db_connection, BioGraphConnection):
            raise TypeError('db_connection must be an instance of the class BioGraphConnection.')
        logging.basicConfig(filename='BiomeDB.log',
                            level=logger_level,
                            format='%(asctime)s %(message)s - %(module)s.%(funcName)s',
                            datefmt='%H:%M:%S-%d.%m.%y')
        self._logger = logging.getLogger(__name__)
        self._logger.info('Initialization of Blaster')
        self.db_connection = db_connection
        self.db_link = self.db_connection.db_link
        self._db_nodes = self._get_db_nodes()
        neo4j._add_header('X-Stream', 'true;format=pretty')
        self.seq_node_dict = {}

    def _get_db_nodes(self):
        db_nodes_dict = {}
        for db_name in ['GenBank', 'PDB', 'UniProt', 'GI']:
            db = self.find_nodes('DB', ['name'], [db_name])
            if not db:
                db, = self.db_connection.data_base.create({'name': db_name})
                db.add_labels('DB')
                db_nodes_dict[db_name] = db
            elif len(db) > 1:
                log_message = 'There are two "%s" database nodes.' % (db_name)
                self._logger.error(log_message)
                raise Warning(log_message)
            else:
                db_nodes_dict[db_name] = db[0]
        return db_nodes_dict

    def find_nodes(self, label, keys=[], values=[]):
        if len(keys) != len(values):
            raise UserWarning('The number of keys must be the same as the number of values.')

        session = cypher.Session(self.db_link)
        transaction = session.create_transaction()
        if keys:
            query = 'MATCH (node:%s) WHERE ' % label
            for key, value in zip(keys, values):
                query += 'node.%s="%s" AND ' % (key, value)
            query = query[:-4] + 'RETURN node'
        else:
            query = 'MATCH (node:%s) RETURN node' % label
        transaction.append(query)
        node = transaction.commit()
        try:
            return list(node[0][0])
        except:
            return []

    def find_organisms(self):
        return self.find_nodes('Organism')

    def _find_node_by_id(self, poly_id):
        try:
            return [self.db_connection.data_base.node(poly_id)]
        except:
            log_message = 'There is no node with id: %s' % poly_id
            print log_message
            self._logger.error(log_message)
            warnings.warn(log_message)
            return []

    def _find_poly_by_id_and_check_seq(self, poly_id, query_seq, long_name_flag):
        poly = self._find_node_by_id(poly_id)
        if not long_name_flag:
            if poly:
                if not poly[0].get_properties()['seq'] == query_seq[:-1]:
                    log_message = 'The sequence from UBLAST result does not match the sequence from the DB! '
                    print log_message
                    self._logger.error(log_message)
                    warnings.warn(log_message)
                else:
                    return poly[0]
            else:
                return []
        else:
            return poly[0]

    # def _line_distinguisher_usearch(self, line):
    #     poly_id, poly_info, query_org, identity, target_seq, target_ref = line.split('\t')[:6]
    #     query_seq = line.split('|')[-1]
    #     return poly_id, poly_info, identity, target_seq.upper(), target_ref, query_seq.upper()

    # def _line_distinguisher_usearch(self, line):
    #     for ublast version 8, must be fixed
    #     head, identity, target_seq, target_ref = line.split('\t')[:4]
    #     poly_id, poly_info, query_org = head.split('|')
    #     query_seq = line.split('\t')[-1]
    #     return poly_id, poly_info, identity, target_seq.upper(), target_ref, query_seq.upper()

    def _line_distinguisher_usearch(self, line):
        # for ublast version 7
        poly_id, poly_info, identity, target_seq, target_ref = line.split('\t')[:5]
        target_seq = target_seq
        evalue, query_seq = line.split('\t')[-2:]
        query_seq = query_seq[:-1]
        return poly_id, poly_info, identity, target_seq.upper(), target_ref, query_seq.upper(), eval(evalue)

    # def _line_distinguisher_usearch(self, line):
    #     # for ublast version 7
    #     poly_id, poly_info, identity, target_seq, target_ref = line.split('\t')[:5]
    #     query_seq = line.split('\t')[-1]
    #     return poly_id, poly_info, identity, target_seq.upper(), target_ref, query_seq.upper()

    def upload_batch_nodes_usearch_sp(self, usearch_result):
        # Read usearch output file
        res_file = open(usearch_result, 'r')
        file_read = res_file.readlines()
        res_file.close()

        # Create a batch
        batch = neo4j.WriteBatch(self.db_connection.data_base)
        line_counter = 0
        # Read file line by line
        for line in file_read:
            if line_counter%100 == 0:
                print 'Processing line %d' % line_counter
            line_counter += 1
            # Distinguish line
            poly_id, poly_info, identity, target_seq, target_ref, query_seq, evalue = self._line_distinguisher_usearch(line)

            # # Check the length of the query sequence
            # if len(query_seq) < 4096:
            #     long_name_flag = False
            # else:
            #     long_name_flag = True
            #     log_message = 'Too long string! String number:%d' % (file_read.index(line))
            #     self._logger.error(log_message)
            #     warnings.warn(log_message)

            # Find poly and check its sequence
            poly = self._find_poly_by_id_and_check_seq(poly_id, query_seq, long_name_flag)

            # If there is no poly, log an error
            if not poly:
                g_start, g_end, ccp = poly_info.split(':')
                log_message = 'Nothing was found by id:%s or sequence:%s ' \
                              'additional info gene_start:%s, gene_end:%s, in ccp:%s' %\
                              (poly_id, query_seq, g_start, g_end, ccp)
                print log_message
                warnings.warn(log_message)
                self._logger.error(log_message)
            else:
                # check if this poly already has an xref
                transaction_out = self._check_xref(target_ref)
                if not transaction_out:
                    # If there is no such pattern, create it: poly--xref--db
                    self._write2batch(batch, target_ref, target_seq, poly, identity, evalue)
                else:
                    # If there is a pattern with a b_poly get its sequence
                    b_poly = transaction_out[0][0]
                    if not long_name_flag:
                        poly_seq = b_poly.get_properties()['seq']

                        # Compare poly's sequence and the file sequence
                        if not poly_seq == target_seq:
                            # If does not match log an error
                            log_message = 'Sequence of the found polypeptide does not match to the one existing in the data base ' \
                                          'ref:%s, seq:%s' % (target_ref, target_seq)
                            print log_message
                            self._logger.error(log_message)

                    # Check if there is already edge 'SIMILAR' between poly and b_poly
                    rels = list(self.db_connection.data_base.match(start_node=b_poly, end_node=poly))
                    if not rels:
                        # If there is no edge 'SIMILAR' between poly and b_poly, create it
                        batch.create(rel(b_poly, ('SIMILAR', {'identity': identity, 'evalue': evalue}), poly))
        batch.submit()

    def _check_xref(self, target_ref):
        # Create cypher session to search pattern
        session = cypher.Session(self.db_link)
        transaction = session.create_transaction()
        query = 'MATCH (p:Polypeptide)--(x:XRef)--(db:DB) ' \
                'WHERE db.name="%s" ' \
                'AND x.id="%s"' \
                'RETURN p' % ('UniProt', target_ref)
        transaction.append(query)
        return transaction.commit()[0]

    def _write2batch(self, batch, target_ref, target_seq, poly, identity, evalue):
        ref = batch.create(node({'id': target_ref}))
        batch.add_labels(ref, 'XRef')
        if not long_name_flag:
            b_poly = batch.create(node({'seq': target_seq}))
        else:
            b_poly = batch.create(node({'seq': ''}))
        batch.add_labels(b_poly, 'Polypeptide', 'BioEntity')
        batch.create(rel(ref, 'LINK_TO', self._db_nodes['UniProt']))
        batch.create(rel(b_poly, 'EVIDENCE', ref))
        batch.create(rel(b_poly, ('SIMILAR', {'identity': identity, 'evalue': evalue}), poly))


class BlastUploader2(BlastUploader):

    def upload_batch_nodes_usearch_sp(self, usearch_result):
        # Read usearch output file
        res_file = open(usearch_result, 'r')
        file_read = res_file.readlines()
        res_file.close()

        # Create a batch
        batch = neo4j.WriteBatch(self.db_connection.data_base)
        line_counter = 0
        # Read file line by line
        for line in file_read:
            if line_counter % 100 == 0:
                print 'Processing line %d' % line_counter
            line_counter += 1
            # Distinguish line
            poly_id, poly_info, identity, target_seq, target_ref, query_seq, evalue, md5_query, md5_target =\
                self._line_distinguisher_usearch(line)

            # Find poly and check its sequence
            poly = self._find_poly_by_id_and_check_seq(poly_id, query_seq)

            # If there is no poly, log an error
            if not poly:
                g_start, g_end, ccp = poly_info.split(':')
                log_message = 'Nothing was found by id:%s or sequence:%s ' \
                              'additional info gene_start:%s, gene_end:%s, in ccp:%s' %\
                              (poly_id, query_seq, g_start, g_end, ccp)
                print log_message
                warnings.warn(log_message)
                self._logger.error(log_message)
            else:
                # try to find Sequence node with query_seq
                try:
                    # try ro find in the dictionary
                    query_seq_node = self.seq_node_dict[md5_query]
                except:
                    # check node with query_md5
                    query_seq_node = list(self.db_connection.data_base.find('Sequence', 'md5', md5_query))
                    if not query_seq_node:
                        # if does not exist, create (seq)-IS_A-(poly)
                        query_seq_node = self._create_seq(batch, query_seq, md5_query)
                        self.seq_node_dict[md5_query] = query_seq_node
                        batch.create(rel(poly, 'IS_A', query_seq_node))
                        batch.set_property(poly, 'md5', md5_query)
                    else:
                        # if exists, do nothing, get node with query_md5
                        query_seq_node = query_seq_node[0]
                        self.seq_node_dict[md5_query] = query_seq_node
                # try to find Sequence node with target_seq
                try:
                    target_seq_node = self.seq_node_dict[md5_target]
                except:
                    # check node with target_md5
                    target_seq_node = list(self.db_connection.data_base.find('Sequence', 'md5', md5_target))
                    if not target_seq_node:
                        # if does not exist, create (seq)-SIMILAR-(seq)-IS_A-(b_poly)--XRef
                        target_seq_node = self._create_seq(batch, target_seq, md5_target)
                        self.seq_node_dict[md5_target] = target_seq_node
                        self._create_similarity_pattern(batch, target_seq_node, query_seq_node, evalue, identity, target_ref, md5_target)
                    else:
                        target_seq_node = target_seq_node[0]
                        self.seq_node_dict[md5_target] = target_seq_node
                batch.create(rel(target_seq_node,
                                 ('SIMILAR', {'identity': identity, 'evalue': evalue}),
                                  query_seq_node))
                    # Possible conditions
                    # if exists,check XRef.id with target_ref
                            # if exists connect node with target_md5: (seq)-SIMILAR-(seq)
                            # if does not exist, connect node with target_md5: (seq)-SIMILAR-(seq)-IS_A-(b_poly)--XRef
        batch.submit()

    def _create_seq(self, batch, seq, md5):
        seq_node = batch.create(node({'md5': md5,
                                      'seq': seq}))
        batch.add_labels(seq_node, 'Sequence')
        return seq_node

    def _create_similarity_pattern(self, batch, target_seq_node, query_seq_node, evalue, identity, target_ref, md5_target):
        b_poly = batch.create(node({'md5': md5_target}))
        batch.add_labels(b_poly, 'Peptide', 'Polypeptide')
        xref = batch.create(node({'id': target_ref}))
        batch.add_labels(xref, 'XRef')
        batch.create(rel(b_poly, 'EVIDENCE', xref))
        batch.create(rel(b_poly, 'IS_A', target_seq_node))
        # batch.create(rel(target_seq_node,
        #                  ('SIMILAR', {'identity': identity, 'evalue': evalue}),
        #                  query_seq_node))


    def _find_poly_by_id_and_check_seq(self, poly_id, query_seq):
        poly = self._find_node_by_id(poly_id)
        if poly:
            if not poly[0].get_properties()['seq'] == query_seq:
                log_message = 'The sequence from UBLAST result does not match the sequence from the DB! '
                print log_message
                self._logger.error(log_message)
                warnings.warn(log_message)
                return []
            else:
                return poly[0]
        else:
            return []

    def _line_distinguisher_usearch(self, line):
        # for ublast version 7
        poly_id, poly_info, identity, target_seq, target_ref = line.split('\t')[:5]
        evalue, query_seq = line.split('\t')[-2:]
        query_seq = query_seq[:-1]
        md5_target = hashlib.md5(target_seq.upper()).hexdigest()
        md5_query = hashlib.md5(query_seq.upper()).hexdigest()
        return poly_id, poly_info, identity, target_seq.upper(), target_ref, query_seq.upper(), eval(evalue), md5_query, md5_target

    def _check_xref(self, target_ref, md5):
        # Create cypher session to search pattern
        session = cypher.Session(self.db_link)
        transaction = session.create_transaction()
        query = 'MATCH (s:Sequence)--(x:XRef)--(db:DB) ' \
                'WHERE db.name="%s" ' \
                'AND x.id="%s" ' \
                'AND s.md5=%.2e' \
                'RETURN s' % ('UniProt', target_ref, md5)
        transaction.append(query)
        return transaction.commit()[0]

    def _check_sequence(self, poly, md5):

        # len(list(graph_db.match(start_node=start_node, end_node=end_node, rel_type=relationship))) > 0:
        # Create cypher session to search pattern
        session = cypher.Session(self.db_link)
        transaction = session.create_transaction()
        query = 'MATCH (s:Sequence)--(x:XRef)--(db:DB) ' \
                'WHERE db.name="%s" ' \
                'AND x.id="%s" ' \
                'AND s.md5=%.2e' \
                'RETURN s' % ('UniProt', target_ref, md5)
        transaction.append(query)
        return transaction.commit()[0]

    def _write2batch(self, batch, target_ref, target_seq, poly, identity, evalue, md5):
        ref = batch.create(node({'id': target_ref}))
        batch.add_labels(ref, 'XRef')
        seq_node = batch.create(node({'seq': target_seq,
                                      'md5': md5,
                                      'evalue': evalue}))
        batch.add_labels(seq_node, 'Sequence', 'AA_Sequence')
        b_poly = batch.create(node({'md5': md5}))
        batch.add_labels(b_poly, 'Polypeptide', 'BioEntity')
        batch.create(rel(ref, 'LINK_TO', self._db_nodes['UniProt']))
        batch.create(rel(seq_node, 'EVIDENCE', ref))
        batch.create(rel(b_poly, ('SIMILAR', {'identity': identity}), poly))
        batch.create(rel(b_poly, 'IS_A', seq_node))
        batch.create(rel(ref, 'LINK_TO', self._db_nodes['UniProt']))

    def _write2batch_add_xref(self, batch, target_ref, md5, identity, poly, seq_node):
        b_poly = self.find_nodes('Polypeptide', ['md5'], md5)
        ref = batch.create(node({'id': target_ref}))
        batch.add_labels(ref, 'XRef')
        batch.create(rel(seq_node, 'EVIDENCE', ref))
        batch.create(rel(b_poly, ('SIMILAR', {'identity': identity}), poly))
        batch.create(rel(ref, 'LINK_TO', self._db_nodes['UniProt']))