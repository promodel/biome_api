import genbank as gb

# gb_file = '/home/artem/work/Biolum/biolum_bacs/full_genomes/Aliivibrio salmonicida LFI1238 plasmid pVSAL840, complete sequence.gb'
# gb_file = 'Chlamydia_trachomatis_A2497_complete_genome_test.gb'
# gb_file = 'Chlamydia_trachomatis_A2497_complete_genome_ver1.gb'
# gb_file = '/home/artem/work/reps/GenBank/Aliivibrio salmonicida LFI1238 chromosome 2.gb'
files = ['Bacillus subtilis subsp. subtilis str. 168 complete genome.gb',
         'Chlamydia trachomatis L2bUCH-1proctitis chromosome, complete genome.gb',
         'Corynebacterium glutamicum ATCC 13032, IS fingerprint type 4-5, complete genome.gb',
         'Escherichia coli str. K-12 substr. MG1655, complete genome.gb']
connection = gb.BioGraphConnection('http://localhost:9494/db/data/')
for gb_file in files:
    print 'Uploading %s' % gb_file
    organism = gb.GenBank(gb_file, connection)
    organism.upload()
print 'All organism has been uploaded.'