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