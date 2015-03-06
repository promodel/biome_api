import os
import genbank as gb
os.chdir = ('/home/artem/work/reps/GenBank/biome_api/biome/load/genbank/')
gb_file_shw = 'Shewanella_woodyi_ATCC_51908_test.gb'
gb_file_chl = 'Chlamydia_trachomatis_A2497_complete_genome_test.gb'
gb_file_ali = 'Aliivibrio_salmonicida_LFI1238_chromosome_2.gb'

connection = gb.BioGraphConnection()
organism = gb.GenBank(gb_file_shw, connection)

gene = list(organism.db_connection.data_base.find('Gene'))[2]
xrefs = organism.rec.features[5].qualifiers['db_xref']
organism.create_xref(xrefs, gene, True)