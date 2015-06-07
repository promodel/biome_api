def parse_mitab(mitab_file):
    open_mitab_file = open(mitab_file, 'r')
    text = open_mitab_file.readlines()
    open_mitab_file.close()
    headers = text[0].split('\n')
    values = text[1].split('\n')
    table_of_headers = headers[0].split('\t')
    table_of_values = values[0].split('\t')
    for header, value in zip(table_of_headers, table_of_values):
        print header, '\n', value, '\n'
#         to be implemented
#         compare with big dictionary and upload current data to DB

def get_mi_dict(dict_file):
    open_file = open(dict_file, 'r')
    text = open_file.readlines()
    open_file.close()
    text = text[text.index('[Term]\n'):]
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
        #     to be implemented
        #     create big dictionary
        text = text[end_index:]