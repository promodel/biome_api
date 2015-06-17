from genbank_ublaster import *
import os

import xml.dom.minidom
from xml.dom.minidom import parse

def parse_mitab(mitab_file):
    # read mitab file
    open_mitab_file = open(mitab_file, 'r')
    text = open_mitab_file.readlines()
    open_mitab_file.close()
    # split headers of the table and data
    headers = text[0].split('\n')
    values = text[1].split('\n')
    table_of_headers = headers[0].split('\t')
    table_of_values = values[0].split('\t')
    # parse and upload it
    for header, value in zip(table_of_headers, table_of_values):
        print header, '\n', value, '\n'
#         to be implemented
#         compare with big dictionary and upload current data to DB

def get_mi_dict(dict_file):
    # read format
    open_file = open(dict_file, 'r')
    text = open_file.readlines()
    open_file.close()
    # cut the description
    text = text[text.index('[Term]\n'):]
    interaction_dict = {}
    # read short parts for each term
    while text[0] == '[Term]\n':
        end_index = text.index('\n')+1
        part = text[0:end_index]
        part_dict = {}
        for line in part:
            try:
                key, value = line.split(': ')
                value = value[:-1]
                if part_dict.has_key(key):
                    part_dict[key] = [part_dict[key], value]
                else:
                    part_dict[key] = value
            except:
                continue
            interaction_dict[part_dict['id']] = part_dict
        #     to be implemented
        #     create big dictionary
        text = text[end_index:]
    return interaction_dict

def filter_eukaryotic(intact_file, substrings=['human', 'eukaryotic'], filtered_filename='intact_filtered.txt'):
    open_file = open(intact_file, 'r')
    filtered_file = open(filtered_filename, 'w')
    for line in open_file:
        if not any(substring in line for substring in substrings):
            filtered_file.write(line)
    open_file.close()
    filtered_file.close()

def cut_mitab(intact_file, short_intact='intact_short.txt', string_number=10):
    open_file = open(intact_file, 'r')
    short_file = open(short_intact, 'w')
    for i in xrange(string_number):
        short_file.write(open_file.next())
    open_file.close()
    short_file.close()

def get_taxons(intact_file):
    open_file = open(intact_file, 'r')
    open_file.readline()
    intact_taxids = {line.split('\t')[28].split(':')[1].split('(')[0] for line in open_file}
    open_file.close()
    return intact_taxids

def compare_taxons(intact_files_directory=['/home/artem/BLAST_DB/intact/pmidMITAB27/'],
                   connection='http://217.26.19.154:17474/db/data/'):
    bgc=BioGraphConnection(connection)
    taxids = set()
    compare_dict = {}
    bu = BlastUploader(bgc, logging.CRITICAL)
    # counter = 0
    for directory in intact_files_directory:
        files = os.listdir(directory)
        try:
            files.remove('.directory')
        except:
            pass
        print 'Directory %s contains %d files' % (directory, len(files))
        for f in files:
            taxids = taxids.union(get_taxons(directory + f))
            # counter += 1
            # if counter%100 == 0:
            #     print len(taxids)
        print 'All files have been read.'
    for taxid in taxids:
        res = bu.find_nodes('Taxon', ['tax_id'], [taxid])
        if res:
            compare_dict[taxid] = 1
        else:
            compare_dict[taxid] = 0
    return compare_dict, taxids

class PSIMIXML():
    def __init__(self):
        self.xml_file = ''
        self. entries_list = ['source', 'experimentList', 'interactorList', 'interactionList']

    def read(self, xml_file='/home/artem/BLAST_DB/intact/1173018.xml'):
        self.xml_file = xml_file
        tree = parse(self.xml_file)
        collection = tree.documentElement
        reactions = collection.getElementsByTagName('entry')
        for reaction in reactions:
            sources = reaction.getElementsByTagName(self.entries_list[0])
            experiments = reaction.getElementsByTagName(self.entries_list[1])
            interactors = reaction.getElementsByTagName(self.entries_list[2])
            interactions = reaction.getElementsByTagName(self.entries_list[3])
        return sources[0], experiments[0], interactors[0], interactions[0]

    def read_interactors(self, interactors):
        interactors_dict = {}
        for interactor in interactors.getElementsByTagName('interactor'):
            # make dictionary of interactors with ids as he keys
            current_interactor = interactor.getAttribute('id')
            interactors_dict[current_interactor] = {}
            # add xrefs to this dictionary
            xrefs = {}
            xrefs[interactor.getElementsByTagName('primaryRef')[0].getAttribute('db')] =\
                interactor.getElementsByTagName('primaryRef')[0].getAttribute('id')
            # for secondary_xref in interactor.getElementsByTagName('secondaryRef'):
            #     xrefs[secondary_xref.getAttribute('db')] =\
            #         secondary_xref.getAttribute('id')
            interactors_dict[current_interactor]['xref'] = xrefs

            # add organism ncbi id
            interactors_dict[current_interactor]['org_id'] = interactor.getElementsByTagName('organism')[0].getAttribute('ncbiTaxId')

            # add sequence
            interactors_dict[current_interactor]['sequence'] = interactor.getElementsByTagName('sequence')[0].firstChild.data

            # add interactor id for crossreferecing in th xml file
            interactors_dict[current_interactor]['interactor_id'] = interactor.getAttribute('id')
        return interactors_dict

    def read_interactions(self, interactions):
        interactions_dict = {}
        for interaction in interactions.getElementsByTagName('interaction'):
            # make dictionary of interactors with ids as he keys
            # current_interaction = interaction.getAttribute('id')
            current_interaction = interaction.getAttribute('imexId')
            interactions_dict[current_interaction] = {}
            # add xrefs to this dictionary
            xrefs = {}
            xrefs[interaction.getElementsByTagName('primaryRef')[0].getAttribute('db')] =\
                interaction.getElementsByTagName('primaryRef')[0].getAttribute('id')
            interactions_dict[current_interaction]['xref'] = xrefs
            # create dictionary for reaction interactor
            participants = interaction.getElementsByTagName('participant')
            interactions_dict[current_interaction]['participants'] = []
            # create experimental ref for xml cross referencing
            experiment_ref = interaction.getElementsByTagName('experimentRef')[0].firstChild.data
            interactions_dict[current_interaction]['exp_ref'] = experiment_ref
            for participant in participants:
                # add xml-cross reference for interactor
                interactions_dict[current_interaction]['participants'].append(
                    participant.getElementsByTagName('interactorRef')[0].firstChild.data)
                # add biological role
                interactions_dict[current_interaction]['bio_role'] = {}
                interactions_dict[current_interaction]['bio_role'] =\
                    participant.getElementsByTagName('biologicalRole')[0].getElementsByTagName('fullName')[0].firstChild.data
        return interactions_dict

    def read_experiments(self, experiments):
        experiments_dict = {}
        for experiment in experiments.getElementsByTagName('experimentDescription'):
            # get experiment id for xml cross referencing
            current_experiment = experiment.getAttribute('id')
            experiments_dict[current_experiment] = {}
            # get xrefs
            xrefs = {}
            xref = experiment.getElementsByTagName('xref')[0]
            xrefs[xref.getElementsByTagName('primaryRef')[0].getAttribute('db')] \
                = xref.getElementsByTagName('primaryRef')[0].getAttribute('id')
            experiments_dict[current_experiment]['xref'] = xrefs
            # get host organism
            experiments_dict[current_experiment]['host_org'] =\
                experiment.getElementsByTagName('hostOrganism')[0].getAttribute('ncbiTaxId')
            #  get interaction detection method
            experiments_dict[current_experiment]['int_detection'] =\
                experiment.getElementsByTagName('interactionDetectionMethod')[0].getElementsByTagName('fullName')[0].firstChild.data
            # get participant identification method
            experiments_dict[current_experiment]['part_identification'] =\
                experiment.getElementsByTagName('participantIdentificationMethod')[0].getElementsByTagName('fullName')[0].firstChild.data
        return experiments_dict

# DOMTree = xml.dom.minidom.parse('/home/artem/BLAST_DB/intact/ecoli/ecoli_01.xml')
# collection=DOMTree.documentElement
# s = set()
# for e in collection.getElementsByTagName('hostOrganism'):
#     s.add(e.getAttribute('ncbiTaxId'))
# ce1=e1.childNodes[0]