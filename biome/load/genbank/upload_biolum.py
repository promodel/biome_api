import genbank as gb
import os
import time

start_time = time.time()

connection = gb.BioGraphConnection()
dir = '/home/artem/work/reps/biolum2/'
files = sorted(os.listdir(dir))

for file in files:
     if '.' in file[0]:
         files.pop(files.index(file))

for i, file in enumerate(files):
    start_file_time = time.time()
    organism = gb.GenBank(dir + file, connection)
    organism.upload()
    print i, file, time.time() - start_file_time
print time.time() - start_time