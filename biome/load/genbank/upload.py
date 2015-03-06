import genbank as gb

gb_file = '/home/artem/work/Biolum/biolum_bacs/full_genomes/Aliivibrio salmonicida LFI1238 plasmid pVSAL840, complete sequence.gb'
connection = gb.BioGraphConnection()
organism = gb.GenBank(gb_file, connection)
organism.upload()