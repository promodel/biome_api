#!/usr/bin/bash
blastp -db $2 -evalue 0.000010 -out "Chlamydia trachomatis A2497_input_blast_blast_out_part$1.xml" -query "Chlamydia trachomatis A2497_input_blast_part$1.FASTA" -outfmt 5 