from ...api import *
from ...load import *
import warnings
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
        if set(list(new_reference) + list(old_reference)) >= bases:
            warnings.warn('Sequences have noncanonical nucleotides!')

        self.new_ref = new_reference
        self.old_ref = old_reference
        self.path = path

    def __repr__(self):
        return "An object of NewReference class"

    def __str__(self):
        return "An object of NewReference class"

    def update_location(self, start, end, flank=4, threshold=8):
        #start -= 1
        if start < 0:
            raise ValueError('Start position must be a positive integer!')
        if end > len(self.old_ref):
            raise ValueError('Site coordinates are out of sequence length!')
        seq = self.old_ref[(start-1):end]

        if len(seq) < threshold:
            # if subsequence is too small
            return self.identify_smallseq_location(start,
                                                   end, flank=flank)
        else:
            return self.identify_seq_location(seq)

    def identify_seq_location(self, seq):
        new_starts = [m.start() + 1 for m in re.finditer('(?=%s)'
                                                         % seq, self.new_ref)]
        if not new_starts:
            return [], [], False
        else:
            return new_starts, [(x + len(seq) - 1) for x in new_starts], True

    def identify_smallseq_location(self, start, end, flank):
        doubleflank = 2*flank
        if (start - flank) >= 0 and (end + flank) <= len(self.old_ref):
            new_start = start - flank
            new_end = end + flank
            seq = self.old_ref[new_start:new_end]
            location = self.identify_seq_location(seq)
            if location[2]:
                return [(x + flank - 1) for x in location[0]], \
                       [x - flank for x in location[1]], True
            else:
                return location
        elif start - doubleflank >= 0:
            new_start = start - doubleflank
            seq = self.old_ref[new_start:end]
            location = self.identify_seq_location(seq)
            if location[2]:
                return [(x + doubleflank - 1) for x in location[0]], \
                       location[1], True
            else:
                return location
        elif end + doubleflank <= len(self.old_ref):
            new_end = end + doubleflank
            seq = self.old_ref[start:new_end]
            location = self.identify_seq_location(seq)
            if location[2]:
                return [x - 1 for x in location[0]], \
                       [x - doubleflank for x in location[1]],\
                       True
            else:
                return location
        else:
            raise ValueError('The threshold is too high!')


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
            location = self.update_location(node.start, node.end, 7, 10)
            if location[2]:
                if len(location[0]) == 1:
                    # easy case! updating the node
                    node.start, node.end = location[0][0], location[1][0]
                    updated += 1
                else:
                    # looking for the closest match
                    diff = [abs(l - node.start) for l in location[0]]
                    i = diff.index(min(diff))
                    pos_loc = ', '.join(['[%d, %d]'
                                         % (location[0][n], location[1][n])
                                         for n in xrange(0, len(location[0]))])

                    logging.warning('Updating problem: there are a few '
                                    'possible locations for %s!\n'
                                    'They are: %s. \n'
                                    'The closest match is [%d, %d], '
                                    'difference = %d.'
                                    % (node, pos_loc,
                                       location[0][i],
                                       location[1][i],
                                       diff[i]))
                    # updating the node
                    node.start, node.end = location[0][i], location[1][i]
                    node.labels += ':ToCheck'
                    node.checking_note = 'A few locations are possible, ' \
                                         'they are %s' % pos_loc
                    several += 1
            else:
                logging.warning('Updating problem: no sequence match'
                                ' for %s!' % node)
                # updating the node
                node.labels += ':ToCheck'
                node.checking_note = 'An old version of the location!'
                notfound += 1

        logging.info('Update of the MetaCyc object for %s was completed!\n'
                     'There were %d updated nodes, %d unchanged nodes, '
                     '%d nodes with a few possible locations.'
                     % (metacyc.organism.name, updated, notfound, several))