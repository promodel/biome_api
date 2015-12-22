# from BiomeDB_classes import *
import logging
import os
from Bio.Blast import NCBIWWW, NCBIXML
from time import time
from py2neo import neo4j, node, rel, cypher
import warnings
import re

class MakeJob():
    def __init__(self, db_connection, e_value=0.00001, logger_level=logging.INFO):
        logging.getLogger("py2neo").setLevel(logging.CRITICAL)
        logging.basicConfig(filename = 'BiomeDB.log',
                            level = logger_level,
                            format = '%(asctime)s %(message)s - %(module)s.%(funcName)s',
                            datefmt='%H:%M:%S-%d.%m.%y')
        self._logger = logging.getLogger(__name__)
        self._logger.info('Initialization of Blaster')
        self._blast_input_txt = ''
        self.e_value = e_value
        self.filename = ''
        self.db_connection = db_connection
        self.db_link = db_connection.db_link

    def split_txt(self, split_number, organism = 'Escherichia_coli_str._K-12_substr._MG1655'):
        #Data base is need to find non analyzed proteins for certain organism.
        if not isinstance(organism, str):
            err_message = 'Organism argument must be a string.'
            self._logger.error(err_message)
            raise ValueError(err_message)

        #Set the name of the file with non-BLASTed  proteins.
        self._blast_input_txt = organism + '_input_blast'
        self._blast_input_txt = re.sub('[!@#$/\s]', '', self._blast_input_txt)
        self.filename = re.sub('[!@#$/]', '', self.filename)
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
                        # for search without taxon
                        header, sequence, tail = lines[line_ind].split('\t')
                        # for search with taxon
                        #header, sequence, tail, org = lines[line_ind].split('|')
                        # org = re.sub('[!@#$/ ]', '_', org)
                        # for search without taxon
                        new_file.write('>' + header + '\t' + tail + '\t' + sequence + '\n')
                        # for search with taxon
                        # new_file.write('>' + header + '|' + tail + '|' + org + sequence + '\n')
                    line_number = average_quantity
                    average_quantity += init_average_quantity
                except:
                    new_file.close()

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
            # for search with taxon
            # transaction_out = self._find_non_blasted_proteins_with_taxon(organism)
            # for search without taxon
            # transaction_out = self._find_non_blasted_proteins(organism)
            # for search for sequences
            transaction_out = self._find_non_blasted_sequences(organism)
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
                        # for search with taxon
                        # poly_str = '%s|%s|%d:%d:%s|%s\n' % (node2link(result[0]), result[1], result[2], result[3], node2link(result[4]), result[5])
                        # for search without taxon
                        poly_str = '%s\t%s\t%d:%d:%s\n' % (node2link(result[0]), result[1], result[2], result[3], node2link(result[4]))
                        polypeptides_file.write(poly_str)
                        poly_counter += 1
                    except:
                        log_message = 'Could not write a sequence into a file: %s' % str(result[0])
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

    def _find_non_blasted_proteins(self, organism):
        session = cypher.Session(self.db_link)
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
        return transaction.commit()[0]

    def _find_non_blasted_sequences(self, organism):
        session = cypher.Session(self.db_link)
        transaction = session.create_transaction()
        query = 'MATCH (s:AA_Sequence)<-[:IS_A]-(p:Polypeptide)<-[:ENCODES]-(g:Gene)-[:PART_OF]->(ccp), (o:Organism) ' \
                'WHERE (p)-[:PART_OF]->(o) ' \
                'AND NOT((s)-[:SIMILAR]-(:AA_Sequence)) ' \
                'AND o.name="%s" ' \
                'AND (ccp)-[:PART_OF]->(o) ' \
                'RETURN distinct s, s.seq, g.start, g.end, ccp' % organism
        transaction.append(query)
        log_message = 'Searching query: ' + query
        self._logger.info(log_message)
        return transaction.commit()[0]

    def _find_non_blasted_proteins_with_taxon(self, organism):
        session = cypher.Session(self.db_link)
        transaction = session.create_transaction()
        query = 'MATCH (t:Taxon)<-[:IS_A]-(o:Organism)<-[:PART_OF]-(p:Polypeptide)<-[:ENCODES]-(g:Gene)-[:PART_OF]->(ccp) ' \
                'WHERE NOT((p)-[:SIMILAR]-(:Polypeptide)) ' \
                'AND o.name="%s" ' \
                'AND HAS(p.seq) AND ' \
                '(ccp)-[:PART_OF]->(o) ' \
                'RETURN distinct p, p.seq, g.start, g.end, ccp, t.scientific_name' % organism
        transaction.append(query)
        log_message = 'Searching query: ' + query
        self._logger.info(log_message)
        return transaction.commit()[0]

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

class MakeJob2(MakeJob):
    def _find_non_blasted_proteins(self, organism):
        session = cypher.Session(self.db_link)
        transaction = session.create_transaction()
        query = 'MATCH (o:Organism)<-[:PART_OF]-(p:Polypeptide)<-[:ENCODES]-(g:Gene)-[:PART_OF]->(ccp) ' \
                'WHERE NOT((p)-[:IS_A]->(:Sequence)) ' \
                'AND o.name="%s" ' \
                'AND HAS(p.seq) AND ' \
                '(ccp)-[:PART_OF]->(o) ' \
                'RETURN distinct p, p.seq, g.start, g.end, ccp' % organism
        transaction.append(query)
        log_message = 'Searching query: ' + query
        self._logger.info(log_message)
        return transaction.commit()[0]