from Bio import SeqIO, SeqRecord

def parse2read(gb_file):
    multi_file = SeqIO.parse(gb_file, 'genbank')
    for i, rec in enumerate(multi_file):
        out_file = open(gb_file + '_part_' + str(i) + '.gb', 'w')
        SeqIO.write(rec, out_file, 'gb')
        out_file.close()