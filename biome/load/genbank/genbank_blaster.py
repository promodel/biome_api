# from BiomeDB_classes import *
import logging
import os
from Bio.Blast import NCBIWWW, NCBIXML
from time import time
from py2neo import neo4j, node, rel, cypher
import warnings

class MakeJob():
    def __init__(self, db_link = 'http://localhost:8484/db/data/', e_value=0.00001, logger_level=logging.INFO):
        logging.basicConfig(filename = 'BiomeDB.log',
                            level = logger_level,
                            format = '%(asctime)s %(message)s - %(module)s.%(funcName)s',
                            datefmt='%H:%M:%S-%d.%m.%y')
        self._logger = logging.getLogger(__name__)
        self._logger.info('Initialization of Blaster')
        self._blast_input_txt = ''
        self.e_value = e_value
        self.filename = ''
        self.db_link = db_link
        self.data_base = neo4j.GraphDatabaseService(self.db_link)

    def split_txt(self, split_number, organism = 'Escherichia_coli_str._K-12_substr._MG1655'):
        #Data base is need to find non analyzed proteins for certain organism.
        if not isinstance(organism, str):
            err_message = 'Organism argument must be a string.'
            self._logger.error(err_message)
            raise ValueError(err_message)

        #Set the name of the file with non-BLASTed  proteins.
        self._blast_input_txt = organism + '_input_blast'
        self.filename = self._blast_input_txt + '.txt'

        #Set current organism
        #Find non-BLASTed proteins or check if the file already exists.
        writing_flag = self._find_non_blasted(organism)

        if writing_flag:
            #Read the file
            poly_file = open(self.filename, 'r')
            lines = poly_file.readlines()
            init_average_quantity = len(lines)/split_number
            average_quantity = init_average_quantity
            line_number = 0

            #Make FASTA-files out of the file
            for i in xrange(split_number):
                new_file = open(self._blast_input_txt + '_part' + str(i) + '.FASTA', 'w')
                if i == split_number - 1:
                    average_quantity = len(lines) + len(lines)%split_number
                try:
                    for line_ind in xrange(line_number, average_quantity):
                            #print line_ind
                            header, sequence, tail = lines[line_ind].split('\t')
                            new_file.write('>' + header + '\t' + tail + sequence + '\n')
                    line_number = average_quantity
                    average_quantity += init_average_quantity
                except:
                    new_file.close()

            # Creating bash-script which will start blastp command
            bash = open('run_blast.sh', 'w')
            bash.write('#!/usr/bin/bash\n')
            self._blast_input_txt + '_part' + str(i) + '.FASTA'
            bash_input = ('blastp -db $2'
                       ' -evalue %f'
                       ' -out "%s_blast_out_part$1.xml"'
                       ' -query "%s' + '_part$1.FASTA"'
                       ' -outfmt 5 '
                       ' -a 4 ') % (self.e_value, self._blast_input_txt, self._blast_input_txt)
            bash.write(bash_input)
            bash.close()

    def _find_non_blasted(self, organism):
        """
        Private method that searches in the data bases for the proteins
        that must be BLASTed. If such file already exs prints message.
        """
        try:
            session = cypher.Session(self.db_link)
        except:
            log_message = 'Could not connect to the base: %s' % self.db_link
            self._logger.error(log_message)
            raise UserWarning(log_message)
        transaction = session.create_transaction()
        query = 'MATCH (org:Organism) WHERE org.name="%s" RETURN org' % organism
        transaction.append(query)
        find_organism = list(transaction.commit())[0]

        if not find_organism:
            err_message = 'There is no organism with name "%s" in the base.' % organism
            self._logger.error(err_message)
            raise ValueError(err_message)
        elif len(find_organism) > 1:
            err_message = 'There are duplicate organisms with name "%s" in the base.' % organism
            self._logger.error(err_message)
            raise UserWarning(err_message)

        if not os.path.isfile(self._blast_input_txt + '.txt'):

            #Find non-BLASTed polypeptides polypeptides for current organism
            transaction = session.create_transaction()
            query = 'MATCH (o:Organism)<-[:PART_OF]-(p:Polypeptide)<-[:ENCODES]-(g:Gene)-[:PART_OF]->(ccp) ' \
                    'WHERE NOT((p)-[:SIMILAR]-(:Polypeptide)) ' \
                    'AND o.name="%s" ' \
                    'AND HAS(p.seq) AND ' \
                    '(ccp)-[:PART_OF]->(o) ' \
                    'RETURN distinct p, p.seq, g.start, g.end, ccp' % organism
            transaction.append(query)
            log_message = 'Searching query: ' + query
            self._logger.info(log_message)
            transaction_out = transaction.commit()[0]
            if not transaction_out:
                log_message = 'Nothing was found on your query: %s' % query
                self._logger.warning(log_message)
                warnings.warn(log_message)
            else:
                #File to write found proteins
                polypeptides_file = open(self.filename, 'w')

                #Write it into the file
                poly_counter = 0
                for result in transaction_out:
                    try:
                        # Write poly-id, seq, start, end, ccp-id
                        poly_str = '%s\t%s\t%d:%d:%s\n' % (node2link(result[0]), result[1], result[2], result[3], node2link(result[4]))
                        polypeptides_file.write(poly_str)
                        poly_counter += 1
                    except:
                        log_message = 'Could not write a polypeptide into a file: %s' % str(result[0])
                        self._logger.error(log_message)
                        warnings.warn(log_message)
                polypeptides_file.close()
                log_message = '%d proteins to be BLASTed.' % poly_counter
                self._logger.info(log_message)
                print log_message
                return True
        else:
            log_message = 'File containing proteins for BLAST is already exist.'
            self._logger.info(log_message)
            print log_message
            return False


def compile_blast_result(xml_file, fasta_file, file_quantity):
    """
    Method reads obtained xml-files and FASTA-files to write final result into txt-file.
    xml_file and fasta_file must have similar name.
    """
    #Open xml-file for reading
    out = open(xml_file + '_parsed.txt', 'w')
    gi_file = open(fasta_file + '_gi.txt', 'w')
    #Read files
    for i in xrange(file_quantity):
        result_file = open(xml_file + '_part' + str(i) + '.xml', 'r')
        initial_file = open(fasta_file + '_part' + str(i) + '.FASTA', 'r')
        blast_read = NCBIXML.parse(result_file)
        #Read result for each protein.
        for reads in blast_read:
            header = initial_file.readline()[1:]
            seq = initial_file.readline()
            out_str = '@\t%s\t%s' % (header[:-1], seq)
            out.write(out_str)
            #Read each alignment for protein
            for alignment in reads.alignments:
                #Write result if it satisfies 2 conditions
                if alignment.hsps[0].identities/float(alignment.length) >= 0.7 and\
                                alignment.length >= round(0.5*len(seq)):
                    partial_alignments = alignment.title.split('>')
                    for partial_alignment in partial_alignments:
                        out_str = '%.2f\t%s\t%s\n' \
                             % (alignment.hsps[0].identities/float(alignment.length),
                                partial_alignment,
                                alignment.hsps[0].sbjct)
                        out.write(out_str)
                    #write protein gi to a file
                    out_str = '%s\n' % (alignment.hit_id.split('|')[1])
                    gi_file.write(out_str)
        result_file.close()
        initial_file.close()
    out.close()
    gi_file.close()


def node2link(gb_node):
    return str(gb_node).split(' ')[0].split('(')[1]


class BlastUploader():
    def __init__(self, db_link='http://localhost:8484/db/data/', logger_level=logging.INFO):
        logging.basicConfig(filename = 'BiomeDB.log',
                            level = logger_level,
                            format = '%(asctime)s %(message)s - %(module)s.%(funcName)s',
                            datefmt='%H:%M:%S-%d.%m.%y')
        self._logger = logging.getLogger(__name__)
        self._logger.info('Initialization of Blaster')
        # self.gdb = Biome_db(db_link)
        self.db_link = db_link
        self.data_base = neo4j.GraphDatabaseService(self.db_link)
        self._gb_node = self._get_gb_node()
        self._db_nodes = self._get_db_nodes()


    def _get_gb_node(self):
        db = list(self.data_base.find('DB', 'name', 'GenBank'))
        if not db:
            db, = self.data_base.create({'name': 'GenBank'})
            db.add_labels('DB')
            return db
        elif len(db) > 1:
            log_message = 'There are two "GenBank" database nodes.'
            self._logger.error(log_message)
            raise Warning(log_message)
        else:
            return db[0]

    def _get_db_nodes(self):
        db_nodes_dict = {}
        for db_name in ['GenBank', 'PDB', 'UniProt', 'GI']:
            db = self.find_nodes('DB', ['name'], [db_name])
            if not db:
                db, = self.data_base.create({'name': db_name})
                db.add_labels('DB')
                db_nodes_dict[db_name] = db
            elif len(db) > 1:
                log_message = 'There are two "%s" database nodes.' % (db_name)
                self._logger.error(log_message)
                raise Warning(log_message)
            else:
                db_nodes_dict[db_name] = db[0]
        return db_nodes_dict

    def find_nodes(self, label, keys=None, values=None):
        # if not (isinstance(keys, list) and isinstance(values, list))\
        #         or not (isinstance(keys, tuple) and isinstance(values, tuple)):
        #     raise UserWarning('Keys and values must be lists or tuples.')
        if len(keys) != len(values):
            raise UserWarning('The number of keys must be the same as the number of values.')

        session = cypher.Session(self.db_link)
        transaction = session.create_transaction()
        query = 'MATCH (node:%s) WHERE ' % label
        for key, value in zip(keys, values):
            query += 'node.%s="%s" AND ' % (key, value)
        query = query[:-4] + 'RETURN node'
        transaction.append(query)
        node = transaction.commit()
        try:
            return list(node[0][0])
        except:
            return []

    def find_organisms(self):
        return self.data_base.find('Organism')

    def _by_id(self, line, method, current_organism=None):
        if method == 'blast':
            node_id = line.split('\t')[1]
        elif method == 'usearch':
            node_id = line.split('\t')[0]
        try:
            return [self.data_base.node(node_id)]
        except:
            log_message = 'There is no node with id: %s' % node_id
            self._logger.error(log_message)
            warnings.warn(log_message)
            return []


    def _by_seq(self, line, method, current_organism):
        if method == 'blast':
            seq = line.split('\t')[-1][:-1]
        elif method == 'usearch':
            seq = line.split('\t')[-1][:-1].upper()
        session = cypher.Session(self.db_link)
        transaction = session.create_transaction()
        query = 'MATCH (o:Organism)<-[:PART_OF]-(p:Polypeptide) ' \
                'WHERE p.seq="%s" ' \
                'AND o.name="%s" ' \
                'RETURN p' % (seq, current_organism)
        transaction.append(query)
        log_message = 'Searching query: ' + query
        self._logger.info(log_message)
        transaction_out = transaction.commit()[0]
        if not transaction_out:
            log_message = 'Nothing was found on the query: '
            self._logger.error(log_message)
            warnings.warn(log_message + query)
            # print log_message, query
            return []
        else:
            return transaction_out[0]

    def terms2dict(self):
        session = cypher.Session(self.db_link)
        transaction = session.create_transaction()
        query = 'MATCH (t:Term) ' \
                'RETURN t, t.text'
        transaction.append(query)
        transaction_out = transaction.commit()[0]
        terms_dict = {}
        for out in transaction_out:
            terms_dict[out[1]] = out[0]
        return terms_dict

    def upload_batch_nodes(self, parsed_result, organism_name=None, seq_flag=True):
        if seq_flag == True and not isinstance(organism_name, str):
            err_message = 'If you want to search polypeptides by sequence,' \
                           ' you have to set the organism name' \
                           ' as e.g. in the name of FASTA-file.'
            self._logger.error(err_message)
            raise ValueError(err_message)

        # Read parsed file
        res_file = open(parsed_result, 'r')
        file_read = res_file.readlines()
        res_file.close()

        start_time = time()
        at_counter = 0

        # Set dictionaries for nodes
        orgs_dict, gi_dict = {}, {}

        # Set dictionaries for batches
        batch_orgs, batch_gi, batch_terms = {}, {}, {}

        # Set dictionary of terms and get terms from the base
        terms_dict = self.terms2dict()

        gi_list, terms_list, orgs_list = [], [], []

        # Choosing the polypeptide nodes searching by id or by sequence
        if seq_flag:
            find_func = self._by_seq
        else:
            find_func = self._by_id

        # Reading file line by line
        batch = neo4j.WriteBatch(self.data_base)
        for line in file_read:
            if line[0] == '@':
                # Polypeptide is found, then searching for homologs
                at_counter += 1

                # Setting current polypeptide
                poly = find_func(line, organism_name)
                if not poly:
                    split_line = line.split('\t')
                    node_id, seq = split_line[1], split_line[-1]
                    g_start, g_end, ccp = split_line[2].split(':')
                    log_message = 'Nothing was found by id:%s or sequence:%s ' \
                                  'additional info gene_start:%s gene_end:%s in ccp:%s' %\
                                  (node_id, seq, g_start, g_end, ccp)
                    warnings.warn(log_message)
                    self._logger.error(log_message)
                batch_out = batch.submit()
                gi_dict, terms_dict, orgs_dict =\
                    self._batch2dict(gi_dict, terms_dict, orgs_dict, gi_list, terms_list, orgs_list, batch_out)
                gi_list, terms_list, orgs_list = [], [], []
                # Creating batch for current poly
                batch_terms, batch_orgs, batch_gi = {}, {}, {}
                batch = neo4j.WriteBatch(self.data_base)
                # print at_counter, time() - start_time, poly
                # start_time = time()
            #     make check of SIMILAR-rel to avoid making replicas
            # Homolog is found
            elif line[4] == '\t' and poly:
                # Reading the homolog line
                similarity, seq, name, org, gi = self._line_distinguisher_blast(line)
                if not gi in gi_dict:
                    # may find more than one polypeptide
                    # maybe we shall include additional property to distinguish polypeptides?
                    # b_poly_search_db = list(self.gdb.data_base.find('Polypeptide', 'seq', seq))

                    # Searching for b_poly
                    b_poly_search_db = self.find_nodes('Polypeptide', ['seq'], [seq])
                    if not b_poly_search_db:
                        # if not gi in batch_gi:
                        b_poly, batch = self._b_poly_process(gi, seq, name, batch_gi, batch)
                        gi_list.append(gi)

                        # Searching for b_poly's organism
                        if not org in orgs_dict:
                            organism_search_db = list(self.data_base.find(
                                'Organism', 'name', org))
                            if not organism_search_db:
                                # if not org in batch_orgs:
                                organism, batch_orgs, batch = self._org_process(org, batch_orgs, batch)
                                orgs_list.append(org)
                            else:
                                organism = organism_search_db[0]
                                orgs_dict[org] = organism
                        else:
                            organism = orgs_dict[org]

                        # Searching for b_poly's term
                        if not name in terms_dict:
                            term_search_db = list(self.data_base.find('Term', 'text', name))
                            if not term_search_db:
                                # if not name in batch_terms:
                                term, batch = self._term_process(name, batch_terms, batch)
                                terms_list.append(name)
                            else:
                                term = term_search_db[0]
                                terms_dict[name] = term
                        else:
                            term = terms_dict[name]

                        # Creating relations
                        ref = batch.create(node({'id': gi}))
                        batch.add_labels(ref, 'XRef')
                        batch.create(rel(ref, 'LINK_TO', self._gb_node))
                        batch.create(rel(b_poly, 'EVIDENCE', ref))
                        batch.create(rel(b_poly, 'PART_OF', organism))
                        batch.create(rel(b_poly, 'HAS_NAME', term))
                        self._similar_process(gi, similarity, poly, b_poly, batch)
                    else:
                        # print at_counter, 'already exists'
                        for b_poly in b_poly_search_db:
                            gi_dict[gi] = b_poly
                            self._similar_process(gi, similarity, poly, b_poly, batch)
                else:
                    b_poly = gi_dict[gi]
                    self._similar_process(gi, similarity, poly, b_poly, batch)
                    # gi_dict[gi] = b_poly
        batch.submit()
        return gi_dict

    def upload_batch_nodes_blast(self, parsed_blast_result, organism_name=None, seq_flag=True):
        if seq_flag == True and not isinstance(organism_name, str):
            err_message = 'If you want to search polypeptides by sequence,' \
                           ' you have to set the organism name' \
                           ' as e.g. in the name of FASTA-file.'
            self._logger.error(err_message)
            raise ValueError(err_message)

        start_time = time()
        at_counter = 0

        # Method to create variables, read file and define function for identify polypeptides
        orgs_dict, gi_dict,\
        batch_orgs, batch_gi, batch_terms,\
        gi_list, terms_list, orgs_list,\
        find_func, terms_dict, file_read = self._prepare_blast(parsed_blast_result, seq_flag)

        # Reading file line by line
        batch = neo4j.WriteBatch(self.data_base)
        for line in file_read:
            if line[0] == '@':
                # Polypeptide is found, then searching for homologs
                at_counter += 1

                # Setting current polypeptide
                poly = find_func(line, organism_name)
                if not poly:
                    split_line = line.split('\t')
                    node_id, seq = split_line[1], split_line[-1]
                    g_start, g_end, ccp = split_line[2].split(':')
                    log_message = 'Nothing was found by id:%s or sequence:%s ' \
                                  'additional info gene_start:%s gene_end:%s in ccp:%s' %\
                                  (node_id, seq, g_start, g_end, ccp)
                    warnings.warn(log_message)
                    self._logger.error(log_message)
                batch_out = batch.submit()
                gi_dict, terms_dict, orgs_dict =\
                    self._batch2dict(gi_dict, terms_dict, orgs_dict, gi_list, terms_list, orgs_list, batch_out)
                gi_list, terms_list, orgs_list = [], [], []
                # Creating batch for current poly
                batch_terms, batch_orgs, batch_gi = {}, {}, {}
                batch = neo4j.WriteBatch(self.data_base)
                # print at_counter, time() - start_time, poly
                # start_time = time()
            #     make check of SIMILAR-rel to avoid making replicas
            # Homolog is found
            elif line[4] == '\t' and poly:
                # Reading the homolog line
                similarity, seq, name, org, gi = self._line_distinguisher_blast(line)
                if not gi in gi_dict:
                    # may find more than one polypeptide
                    # maybe we shall include additional property to distinguish polypeptides?
                    # b_poly_search_db = list(self.gdb.data_base.find('Polypeptide', 'seq', seq))

                    # Searching for b_poly
                    b_poly_search_db = self.find_nodes('Polypeptide', ['seq'], [seq])
                    if not b_poly_search_db:
                        # if not gi in batch_gi:
                        b_poly, batch = self._b_poly_process(gi, seq, name, batch_gi, batch)
                        gi_list.append(gi)

                        # Searching for b_poly's organism
                        if not org in orgs_dict:
                            organism_search_db = list(self.data_base.find(
                                'Organism', 'name', org))
                            if not organism_search_db:
                                # if not org in batch_orgs:
                                organism, batch_orgs, batch = self._org_process(org, batch_orgs, batch)
                                orgs_list.append(org)
                            else:
                                organism = organism_search_db[0]
                                orgs_dict[org] = organism
                        else:
                            organism = orgs_dict[org]

                        # Searching for b_poly's term
                        if not name in terms_dict:
                            term_search_db = list(self.data_base.find('Term', 'text', name))
                            if not term_search_db:
                                # if not name in batch_terms:
                                term, batch = self._term_process(name, batch_terms, batch)
                                terms_list.append(name)
                            else:
                                term = term_search_db[0]
                                terms_dict[name] = term
                        else:
                            term = terms_dict[name]

                        # Creating relations
                        ref = batch.create(node({'id': gi}))
                        batch.add_labels(ref, 'XRef')
                        batch.create(rel(ref, 'LINK_TO', self._gb_node))
                        batch.create(rel(b_poly, 'EVIDENCE', ref))
                        batch.create(rel(b_poly, 'PART_OF', organism))
                        batch.create(rel(b_poly, 'HAS_NAME', term))
                        self._similar_process(gi, similarity, poly, b_poly, batch)
                    else:
                        # print at_counter, 'already exists'
                        for b_poly in b_poly_search_db:
                            gi_dict[gi] = b_poly
                            self._similar_process(gi, similarity, poly, b_poly, batch)
                else:
                    b_poly = gi_dict[gi]
                    self._similar_process(gi, similarity, poly, b_poly, batch)
                    # gi_dict[gi] = b_poly
        batch.submit()
        return gi_dict

    def _prepare_blast(self, parsed_result, seq_flag):
        # Read parsed file
        res_file = open(parsed_result, 'r')
        file_read = res_file.readlines()
        res_file.close()

        # Set dictionaries for nodes
        orgs_dict, gi_dict = {}, {}

        # Set dictionaries for batches
        batch_orgs, batch_gi, batch_terms = {}, {}, {}

        # Set dictionary of terms and get terms from the base
        terms_dict = self.terms2dict()

        gi_list, terms_list, orgs_list = [], [], []

        # Choosing the polypeptide nodes searching by id or by sequence
        if seq_flag:
            find_func = self._by_seq
        else:
            find_func = self._by_id
        return orgs_dict, gi_dict,\
               batch_orgs, batch_gi, batch_terms,\
               gi_list, terms_list, orgs_list,\
               find_func, terms_dict, file_read

    def _prepare_usearch(self, parsed_result, seq_flag):
        # Read parsed file
        res_file = open(parsed_result, 'r')
        file_read = res_file.readlines()
        res_file.close()

        # Set dictionaries for nodes
        sp_dict = {}

        # Set dictionaries for batches
        batch_sp = {}

        sp_list = []

        # Choosing the polypeptide nodes searching by id or by sequence
        if seq_flag:
            find_func = self._by_seq
        else:
            find_func = self._by_id

        # Get protein's xrefs
        existing_poly_xrefs = self._get_xrefs_for_polypeptides
        return sp_dict, batch_sp, sp_list, find_func, file_read, existing_poly_xrefs

    def upload_batch_nodes_usearch(self, usearch_result, organism_name=None, seq_flag=True):
        ref_dict, batch_sp, sp_list, find_func, file_read, existing_poly_xrefs = self._prepare_usearch(usearch_result, seq_flag)

        at_counter = 0
        batch = neo4j.WriteBatch(self.data_base)
        for line in file_read:
            at_counter += 1

            # Setting current polypeptide
            poly = find_func(line, 'usearch', organism_name)
            poly = poly[0]

            # Reading the usearch line
            poly_id, poly_info, identity, target_seq, target_ref, query_seq = self._line_distinguisher_usearch(line)
            # sp_list, terms_list, orgs_list = [], [], []

            # Creating batch for current poly
            batch = neo4j.WriteBatch(self.data_base)
            if not poly:
                g_start, g_end, ccp = poly_info.split(':')
                log_message = 'Nothing was found by id:%s or sequence:%s ' \
                              'additional info gene_start:%s gene_end:%s in ccp:%s' %\
                              (poly_id, query_seq, g_start, g_end, ccp)
                warnings.warn(log_message)
                self._logger.error(log_message)

            # print poly_info, identity, target_ref
            target_ref, db_name = self._identify_ref(target_ref)
            if not existing_poly_xrefs[db_name].has_key(target_ref):
                # Searching for b_poly
                b_poly_search_db = self.find_nodes('Polypeptide', ['seq'], [target_seq])
                if not b_poly_search_db:
                    # if not gi in batch_gi:

                    # must be a batch-dict for refs
                    b_poly, batch, batch_sp = self._b_poly_process_usearch(target_ref, target_seq, batch, batch_sp)
                    sp_list.append(target_ref)
                    # Creating relations
                    ref = batch.create(node({'id': target_ref}))
                    batch.add_labels(ref, 'XRef')
                    batch.create(rel(ref, 'LINK_TO', self._db_node))
                    batch.create(rel(b_poly, 'EVIDENCE', ref))
                    # sp_dict[target_ref] = b_poly
                    self._similar_process_usearch(target_ref, identity, poly, b_poly, batch)
                else:
                    # print at_counter, 'already exists'
                    for b_poly in b_poly_search_db:
                        ref_dict[target_ref] = b_poly
                        self._similar_process_usearch(target_ref, identity, poly, b_poly, batch)
            else:
                b_poly = ref_dict[target_ref]
                self._similar_process_usearch(target_ref, identity, poly, b_poly, batch)
                # gi_dict[gi] = b_poly
        batch.submit()
        return ref_dict

    def _identify_ref(self, ref):
        signature = ref.split('|')[0]
        db_names = ['GenBank', 'GI', 'UniProt', 'PDB']
        if signature == 'gi':
            db_name = db_names[1]
        elif signature == 'sp':
            db_name = db_names[2]
        self._db_node = self._db_nodes[db_name]
        return ref.split('|')[1], db_name

    def _batch2dict(self, gi_dict, terms_dict, orgs_dict, gi_list, terms_list, orgs_list, batch_res):
        batch_res = [res for res in batch_res if type(res) == neo4j.Node]
        gi_counter, terms_counter, orgs_counter = 0, 0, 0
        for res in batch_res:
            labels = res.get_labels()
            if 'Polypeptide' in labels:
                gi_dict[gi_list[gi_counter]] = res
                gi_counter += 1
            elif 'Organism' in labels:
                orgs_dict[orgs_list[orgs_counter]] = res
                orgs_counter += 1
            elif 'Term' in labels:
                terms_dict[terms_list[terms_counter]] = res
                terms_counter += 1
        return gi_dict, terms_dict, orgs_dict


    def _b_poly_process_usearch(self, target_ref, seq, batch, batch_sp):
        if not batch_sp.has_key(target_ref):
            b_poly = batch.create(node({'seq': seq}))
            batch.add_labels(b_poly, 'Polypeptide', 'BioEntity')
            batch_sp[target_ref] = b_poly
        else:
            b_poly = batch_sp[target_ref]
        return b_poly, batch, batch_sp

    def _b_poly_process(self, gi, seq, name, batch_gi, batch):
        if not gi in batch_gi:
            b_poly = batch.create(node({'seq': seq, 'name': name}))
            batch.add_labels(b_poly, 'Polypeptide', 'BioEntity')
            batch_gi[gi] = b_poly
        else:
            b_poly = batch_gi[gi]
        return b_poly, batch

    def _org_process(self, org, batch_orgs, batch):
        if not org in batch_orgs:
            organism = batch.create(node({'name': org, 'source': ['GenBank']}))
            batch.add_labels(organism, 'Organism')
            batch_orgs[org] = organism
            # batch_orgs.append(organism.body.values()[0])
        else:
            organism = batch_orgs[org]
        return organism, batch_orgs, batch

    def _term_process(self, name, batch_terms, batch):
        if not name in batch_terms:
            term = batch.create(node({'text': name}))
            batch.add_labels(term, 'Term')
            batch_terms[name] = term
            # batch_terms.append(term.body.values()[0])
            # terms_dict[names[org_counter]] = term
        else:
            term = batch_terms[name]
        return term, batch

    def _similar_process_usearch(self, target_ref, identity, poly, b_poly, batch):
        # add e_value
        batch.create(rel(b_poly,
                         ('SIMILAR', {'similarity': identity}),
                         poly))

    def _similar_process(self, gi, similarity, poly, b_poly, batch):
        for multiple_poly in poly:
            batch.create(rel(b_poly,
                             ('SIMILAR', {'similarity': similarity}),
                             multiple_poly))

    def _find_or_create_node(self, key_string, node_dict, batch_dict, node_type, batch):
        if not key_string in node_dict:
            node_search_db = list(self.gdb.data_base.find(node_type, 'name', key_string))
            if not node_search_db:
                if not key_string in batch_dict:
                    current_node = batch.create(node({'name': key_string}))
                    batch.add_labels(current_node, node_type)
                    batch_dict[key_string] = current_node
                    # batch_orgs.append(organism.body.values()[0])
                else:
                    current_node = batch_dict[key_string]
            else:
                current_node = node_search_db[0]
                node_dict[key_string] = current_node

    def _line_distinguisher_blast(self, line):
        similarity, refs, seq = line.split('\t')
        similarity = float(similarity)
        seq = seq[:-1]
        refs = refs[3:]
        refs = refs.split('|')
        gi = refs[0::4][0]
        names_and_orgs = refs[3::4]
        name = [name.split('[')[0][1:-1] for name in names_and_orgs][0]
        org = [org.split('[')[1].split(']')[0] for org in names_and_orgs][0]
        return similarity, seq, name, org, gi

    def _line_distinguisher_usearch(self, line):
        poly_id, poly_info, identity, target_seq, target_ref = line.split('\t')[:5]
        query_seq = line.split('\t')[-1]
        return poly_id, poly_info, identity, target_seq, target_ref, query_seq

    def _get_xrefs_for_polypeptides(self):
        db_xref_list = {}
        session = cypher.Session(self.db_link)
        for db_name in ['GenBank', 'PDB', 'UniProt', 'GI']:
            transaction = session.create_transaction()
            query = 'MATCH (x:XRef)--(db:DB) ' \
                    'WHERE db.name="%s" ' \
                    'RETURN x, x.id' % (db_name)
            transaction.append(query)
            log_message = 'Searching query: ' + query
            self._logger.info(log_message)
            transaction_out = transaction.commit()[0]
            # res[0][0].values[0]
            xref_dict = {}
            for line in transaction_out:
                xref_dict[line[1]] = line[0]
            db_xref_list[db_name] = xref_dict
        return db_xref_list

    def upload_batch_nodes_usearch_sp(self, usearch_result, organism_name=None, seq_flag=True):
        # Read usearch output file
        res_file = open(usearch_result, 'r')
        file_read = res_file.readlines()
        res_file.close()

        # Get uniprot-db node
        sp_node = self._db_nodes['UniProt']

        # Choose the polypeptide nodes searching by id or by sequence
        if seq_flag:
            find_func = self._by_seq
        else:
            find_func = self._by_id

        # Create a batch
        batch = neo4j.WriteBatch(self.data_base)

        # Read file line by line
        for line in file_read:
            # Find poly to attach his homologs
            poly = find_func(line, 'usearch', organism_name)
            poly = poly[0]

            # If there is no poly, log an error
            if not poly:
                g_start, g_end, ccp = poly_info.split(':')
                log_message = 'Nothing was found by id:%s or sequence:%s ' \
                              'additional info gene_start:%s gene_end:%s in ccp:%s' %\
                              (poly_id, query_seq, g_start, g_end, ccp)
                warnings.warn(log_message)
                self._logger.error(log_message)

            # Distinguish line
            poly_id, poly_info, identity, target_seq, target_ref, query_seq = self._line_distinguisher_usearch(line)

            # Create cypher session to search pattern
            session = cypher.Session(self.db_link)
            transaction = session.create_transaction()
            query = 'MATCH (p:Polypeptide)--(x:XRef)--(db:DB) ' \
                    'WHERE db.name="%s" ' \
                    'AND x.id="%s"' \
                    'RETURN p' % ('UniProt', target_ref)
            transaction.append(query)
            transaction_out = transaction.commit()[0]
            if not transaction_out:
                # If there is no such pattern, create it
                ref = batch.create(node({'id': target_ref}))
                batch.add_labels(ref, 'XRef')
                b_poly = batch.create(node({'seq': target_seq}))
                batch.add_labels(b_poly, 'Polypeptide', 'BioEntity')
                batch.create(rel(ref, 'LINK_TO', sp_node))
                batch.create(rel(b_poly, 'EVIDENCE', ref))
                batch.create(rel(b_poly, 'SIMILAR', poly))
            else:
                # If there is a pattern get the sequence of the found poly
                b_poly = transaction_out[0][0]
                poly_seq = b_poly.get_properties()['seq']

                # Compare poly's sequence and the file sequence
                if poly_seq == target_seq:
                    rels = list(self.data_base.match(start_node=b_poly, end_node=poly))
                    if not rels:
                        # If matches, create 'SIMILAR' edge
                        batch.create(rel(b_poly, 'SIMILAR', poly))
                else:
                    # If does not match log an error
                    log_message = 'Sequence if the found polypeptide to the one existing in the data base' \
                                  'ref:%s, seq:%s' % (target_ref, target_seq)
                    self._logger.error(log_message)

        batch.submit()