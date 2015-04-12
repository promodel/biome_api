import subprocess
import os
from ftplib import FTP
import logging
import gzip
import genbank_ublaster
from Bio import SeqIO
import datetime
import sys

class PrepareSwissProt():
    def __init__(self, logger_level=logging.INFO):
        logging.basicConfig(filename='BiomeDB.log',
                            level=logger_level,
                            format='%(asctime)s %(message)s - %(module)s.%(funcName)s',
                            datefmt='%H:%M:%S-%d.%m.%y')
        self._logger = logging.getLogger(__name__)
        self._logger.info('Initialization of PrepareSwissProt')
        self.server_address = 'ftp.uniprot.org'
        self.path = ''
        try:
            self.ftp_connection = FTP(self.server_address)
            self.ftp_connection.login()
            print self.ftp_connection.welcome
        except:
            log_message = 'Could not connect to server %s' % (self.server_address)
            self._logger.error(log_message)
            raise RuntimeError(log_message)

    def change_folder(self, path='/pub/databases/uniprot/current_release/knowledgebase/proteomes'):
        self._refresh_connection()
        self.ftp_connection.cwd(path)
        self.path = path
        print 'Current path:\n%s' % (self.server_address + self.path)

    def list_folder(self):
        self._refresh_connection()
        return self.ftp_connection.nlst()

    def download_file(self, filename, force=False):
        self._refresh_connection()
        if force or not self._check_existing_file(filename) or not self._check_file_date(filename):
            try:
                print '%s to download.' % convert_bytes(self.ftp_connection.size(filename))
                self.ftp_connection.retrbinary('RETR ' + filename, open(filename, 'w').write)
                return True
            except:
                log_message = 'Could not connect to server %s' % (self.server_address)
                self._logger.error(log_message)
                raise UserWarning(log_message)
        else:
            print 'There is the last version of file %s' % filename
            return False

    def _check_existing_file(self, filename):
        if os.path.isfile(filename):
            return True
        else:
            return False

    def _check_file_date(self, filename):
        ftp_time = self.ftp_connection.sendcmd('MDTM ' + filename)
        if ftp_time[:3] == '213':
            local_time = os.path.getatime(filename)
            ftp_time = datetime.datetime.strptime(ftp_time[4:], '%Y%m%d%H%M%S')
            local_time = datetime.datetime.fromtimestamp(local_time)
            if ftp_time > local_time:
                log_message = 'File %s is up to date.' % (filename)
                self._logger.info(log_message)
                return False
            else:
                log_message = 'There is a new version of file %s is up to date.' % (filename)
                self._logger.info(log_message)
                return True
        else:
            log_message = 'Could not connect to server %s' % (self.server_address)
            self._logger.error(log_message)
            raise UserWarning(log_message)

    def _refresh_connection(self):
        try:
            self.ftp_connection.voidcmd('NOOP')
        except:
            self.ftp_connection = FTP(self.server_address)
            self.ftp_connection.login()
            self.ftp_connection.cwd(self.path)
            print 'Connection was refreshed.'

    def unzip_file(self, filename):
        open_archive = gzip.GzipFile(filename, 'rb')
        read_archive = open_archive.read()
        open_archive.close()

        unzipped_file = file(filename[:-3], 'wb')
        unzipped_file.write(read_archive)
        unzipped_file.close()
        print filename + ' was successfully unzipped.'

    def convert_swissprot2fasta(self, swiss_name, fasta_name):
        inputhandle = open(swiss_name,  "rU")# rU = read universal
        outputhandle = open(fasta_name,  "w")# w = write
        try:
            read_swiss = SeqIO.parse(inputhandle,  "swiss")# swiss is for the SwissProt *dat format
            SeqIO.write(read_swiss, outputhandle,  "fasta")
        except:
            inputhandle.close()
            outputhandle.close()
            log_message = 'Could not convert %s' % (swiss_name)
            self._logger.error(log_message)
            raise SystemError(log_message)
        inputhandle.close()
        outputhandle.close()
        print swiss_name + ' was successfully converted to ' + fasta_name
        os.remove(swiss_name + '.dat')

    def close_connection(self):
        self.ftp_connection.close()

    def make_ublast(self, db='/home/artem/BLAST_DB/refseq.udb',
                    input_file='/home/artem/work/reps/GenBank/biome_api/biome/load/genbank/Chlamydia\ trachomatis\ L2bUCH-1proctitis_input_blast_part0.FASTA',
                    output_file='Python_shell_check_usearch_chalmidiya.txt',
                    usearch_location='~/BLAST_DB/usearch7.0.1090_i86linux32'):
        cmd = 'time %s' \
              ' -ublast %s' \
              ' -db %s' \
              ' -evalue 1e-5' \
              ' -mid 70.0' \
              ' -self' \
              ' -minqt 0.5' \
              ' -maxaccepts 500' \
              ' -userout %s' \
              ' -userfields query+id+tseq+target+qs+ts+qlor+qhir+tlor+thir+evalue+qseq' % (usearch_location, input_file, db, output_file)
        process = subprocess.Popen(cmd, shell=True)
        process.wait()
        if not process:
            log_message = 'Could not run UBLAST.'
            self._logger.error(log_message)
            raise SystemError(log_message)
        else:
            print 'UBLAST was successful.'

    def create_udb(self, fastafile, udbfile):
        cmd = '~/BLAST_DB/usearch7.0.1090_i86linux32 -makeudb_ublast %s -output %s' % (fastafile, udbfile)
        process = subprocess.Popen(cmd, shell=True)
        process.wait()
        if not process:
            log_message = 'Could not convert FASTA to UDB.'
            self._logger.error(log_message)
            raise SystemError(log_message)
        else:
            print 'Successfully converted FASTA to UDB.'

def convert_bytes(num):
    for x in ['bytes','KB','MB','GB']:
        if num < 1024.0:
            return "%3.1f %s" % (num, x)
        num /= 1024.0
    return "%3.1f %s" % (num, 'TB')

def main():
    try:
        psp = PrepareSwissProt()
        psp.change_folder('/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/')
        for filename in ['uniprot_sprot_bacteria.dat.gz']:
        # for filename in ['LICENSE']:
            try:
                download_result = psp.download_file(filename, force=False)
                if not download_result:
                    sys.exit(2)
            except:
                print 'File has not been downloaded.'
                sys.exit(2)
            try:
                psp.unzip_file(filename)
            except:
                print 'File has not been unzipped.'
                sys.exit(2)
            try:
                filename = filename.split('.')[0]
                psp.convert_swissprot2fasta(filename + '.dat', filename + '.fasta')
            except:
                print 'File has not been converted from .dat to .fasta'
                sys.exit(2)
            try:
                psp.create_udb(filename + '.fasta', filename + '.udb')
            except:
                print 'File has not been converted from .fasta to .udb'
    except:
        print 'Main function crashed.'
        sys.exit(2)

if __name__ == "__main__":
    main()


# import doctest
# doctest.testfile('usearch_test.txt')