import genbank as gb
import os
import time
from py2neo import cypher


start_time = time.time()
connection = gb.BioGraphConnection('http://localhost:8484/db/data/')
no = gb.GenomeRelations(connection)
orgs = no.find_organisms()
for org in orgs:
    print org
    no.set_organism(int(gb.node2link(org)))
    no.make_next_relation()
    no.make_overlap_relation()
print 'NEXT-OVERLAP', time.time() - start_time

start_time = time.time()
session = cypher.Session(connection.db_link)
transaction = session.create_transaction()
query = 'MATCH (t:Taxon),(o:Organism) WHERE t.scientific_name=o.name ' \
        'CREATE o-[r:IS_A]->t'
transaction.append(query)
transaction.commit()
print 'Taxonomy', time.time() - start_time