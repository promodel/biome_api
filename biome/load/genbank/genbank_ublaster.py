import logging
from time import time
from py2neo import neo4j, node, rel, cypher
import warnings
from genbank import BioGraphConnection
from genbank_blaster import MakeJob

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
        query_seq = line.split('\t')[-1]
        return poly_id, poly_info, identity, target_seq.upper(), target_ref, query_seq.upper()

    # def upload_batch_nodes_usearch_sp(self, usearch_result):
    #     # Read usearch output file
    #     res_file = open(usearch_result, 'r')
    #     file_read = res_file.readlines()
    #     res_file.close()
    #
    #     # Create a batch
    #     batch = neo4j.WriteBatch(self.db_connection.data_base)
    #     line_counter = 0
    #     # Read file line by line
    #     for line in file_read:
    #         if line_counter%100 == 0:
    #             print 'Processing line %d\n' % line_counter
    #         line_counter += 1
    #         # Distinguish line
    #         if not len(line) > 4096:
    #             poly_id, poly_info, identity, target_seq, target_ref, query_seq = self._line_distinguisher_usearch(line)
    #
    #             # Find poly to attach his homologs
    #             poly = self._find_poly_by_id_and_check_seq(poly_id, query_seq)
    #
    #             # If there is no poly, log an error
    #             if not poly:
    #                 g_start, g_end, ccp = poly_info.split(':')
    #                 log_message = 'Nothing was found by id:%s or sequence:%s ' \
    #                               'additional info gene_start:%s, gene_end:%s, in ccp:%s' %\
    #                               (poly_id, query_seq, g_start, g_end, ccp)
    #                 print log_message
    #                 warnings.warn(log_message)
    #                 self._logger.error(log_message)
    #             else:
    #                 # check xref
    #                 transaction_out = self._check_xref(target_ref)
    #                 if not transaction_out:
    #                     # If there is no such pattern, create it
    #                     self._write2batch(batch, target_ref, target_seq, poly)
    #                 else:
    #                     # If there is a pattern get the sequence of the found poly
    #                     b_poly = transaction_out[0][0]
    #                     poly_seq = b_poly.get_properties()['seq']
    #
    #                     # Compare poly's sequence and the file sequence
    #                     if poly_seq == target_seq:
    #                         rels = list(self.db_connection.data_base.match(start_node=b_poly, end_node=poly))
    #                         if not rels:
    #                             # If matches, create 'SIMILAR' edge
    #                             batch.create(rel(b_poly, 'SIMILAR', poly))
    #                     else:
    #                         # If does not match log an error
    #                         log_message = 'Sequence of the found polypeptide does not match to the one existing in the data base ' \
    #                                       'ref:%s, seq:%s' % (target_ref, target_seq)
    #                         print log_message
    #                         self._logger.error(log_message)
    #         else:
    #             poly_id, poly_info, identity, target_seq, target_ref, query_seq = self._line_distinguisher_usearch(line)
    #
    #             log_message = 'Too long string! String number:%d' % (file_read.index(line))
    #             self._logger.error(log_message)
    #             warnings.warn(log_message)
    #     batch.submit()

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
                print 'Processing line %d\n' % line_counter
            line_counter += 1
            # Distinguish line
            poly_id, poly_info, identity, target_seq, target_ref, query_seq = self._line_distinguisher_usearch(line)

            # Check the length of the query sequence
            if len(query_seq) < 4096:
                long_name_flag = False
            else:
                long_name_flag = True
                log_message = 'Too long string! String number:%d' % (file_read.index(line))
                self._logger.error(log_message)
                warnings.warn(log_message)

            # Find poly to attach his homologs
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
                # check xref
                transaction_out = self._check_xref(target_ref)
                if not transaction_out:
                    # If there is no such pattern, create it
                    self._write2batch(batch, target_ref, target_seq, poly, long_name_flag)
                else:
                    # If there is a pattern get the sequence of the found poly
                    b_poly = transaction_out[0][0]
                    if not long_name_flag:
                        poly_seq = b_poly.get_properties()['seq']

                        # Compare poly's sequence and the file sequence
                        if poly_seq == target_seq:
                            rels = list(self.db_connection.data_base.match(start_node=b_poly, end_node=poly))
                            if not rels:
                                # If matches, create 'SIMILAR' edge
                                batch.create(rel(b_poly, 'SIMILAR', poly))
                        else:
                            # If does not match log an error
                            log_message = 'Sequence of the found polypeptide does not match to the one existing in the data base ' \
                                          'ref:%s, seq:%s' % (target_ref, target_seq)
                            print log_message
                            self._logger.error(log_message)
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

    def _write2batch(self, batch, target_ref, target_seq, poly, long_name_flag):
        ref = batch.create(node({'id': target_ref}))
        batch.add_labels(ref, 'XRef')
        if not long_name_flag:
            b_poly = batch.create(node({'seq': target_seq}))
        else:
            b_poly = batch.create(node({'seq': ''}))
        batch.add_labels(b_poly, 'Polypeptide', 'BioEntity')
        batch.create(rel(ref, 'LINK_TO', self._db_nodes['UniProt']))
        batch.create(rel(b_poly, 'EVIDENCE', ref))
        batch.create(rel(b_poly, 'SIMILAR', poly))