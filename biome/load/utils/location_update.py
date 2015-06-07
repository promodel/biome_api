from ...api import *
import warnings
from time import ctime, time
import logging



class NewReference():
    """
    """
    def __init__(self, new_reference, old_reference, path='./'):
        bases = set(['A', 'T', 'G', 'C', 'U', 'W', 'S', 'M', 'K', 'R',
                     'Y', 'B', 'D', 'H', 'V', 'N', '-'])
        if not isinstance(new_reference, basestring):
            raise TypeError('The new_reference argument must be a string!')
        if not isinstance(old_reference, basestring):
            raise TypeError('The old_reference argument must be a string!')
        if not isinstance(path, basestring):
            raise TypeError('The path argument must be a string!')
        if not os.path.isdir(path):
            raise ValueError('The path does not exist!')

        new_reference = new_reference.upper()
        old_reference = old_reference.upper()
        if set(list(new_reference) + list(old_reference)) <= bases:
            warnings.warn('Sequences have noncanonical nucleotides!')

        self.new_ref = new_reference
        self.old_ref = old_reference
        self.path = path

    def __repr__(self):
        return "An object of NewReference class"

    def __str__(self):
        return "An object of NewReference class"

    def identify_location(self, start, end):
        seq = self.old_ref[(start-1):end]
        new_starts = [m.start()+1 for m in re.finditer('(?=%s)'
                                                       % seq, self.new_ref)]
        if not new_starts:
            return [], [], False
        else:
            return new_starts, map(lambda x: x + len(seq) - 1, new_starts), True

    def update_metacyc(self, metacyc, ccp_node):
        if not isinstance(metacyc, MetaCyc):
            raise TypeError("The metacyc argument must be of the "
                            "MetaCyc class!")
        logging.basicConfig(filename='%slocation_update.log' % self.path,
                            level=logging.INFO,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%H:%M:%S-%d.%m.%y')
        logging.info('Updating a MetaCyc object for %s.\n'
                     'Genetic object - %s.'
                     % (metacyc.organism.name, ccp_node.name))

        updated, several, notfound = [0]*3
        edges_ccp = [e.source for e in metacyc.edges
                     if e.target == ccp_node and e.label == 'PART_OF']
        nodes = [n for n in metacyc.genes + metacyc.terminators +
                            metacyc.promoters + metacyc.BSs + metacyc.other_nodes
                 if hasattr(n, 'start') and hasattr(n, 'end') and n in edges_ccp]

        for node in nodes:
            location = self.identify_location(node.start, node.end)
            if location[2] is True:
                if len(location[0]) == 1:
                    node.start, node.end = location[0][0], location[1][0]
                    updated += 1
                else:
                    diff = [abs(l - node.start) for l in location[0]]
                    i = diff.index(min(diff))
                    node.start, node.end = location[0][i], location[1][i]
                    several += 1
                    pos_loc = ', '.join(['[%d, %d]'
                                         % (location[0][n], location[1][n])
                                         for n in xrange(0, len(location[0]))])

                    logging.warning('Updating problem: there are a few '
                                    'possible locations for %s!\n'
                                    'They are: %s. \n'
                                    'The closest match is [%d, %d], '
                                    'difference = %d.'
                                    % (node.__str__, pos_loc,
                                       location[0][i],
                                       location[1][i],
                                       diff[i]))
                    node.labels += ':ToCheck'
                    node.checking_note = 'A few locations are possible, ' \
                                         'they are %s' % pos_loc
            else:
                notfound += 1
                logging.warning('Updating problem: no sequence match'
                                ' for %s!' % node.__str__)
                node.labels += ':ToCheck'
                node.checking_note = 'An old version of the location!'
        logging.info('Update of the MetaCyc object for %s was completed!\n'
                     'There were %d updated nodes, %d unchanged nodes,'
                     '%d nodes with a few possible locations.'
                     % (metacyc.organism.name, updated, notfound, several))




