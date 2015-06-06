from ...api import *
import warnings

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
        new_starts = [m.start()+1 for m in re.finditer('(?=%s)' % seq, self.new_ref)]
        if not new_starts:
            return [], [], False
        else:
            return new_starts, map(lambda x: x + len(seq) - 1, new_starts), True

    def update_metacyc(self, metacyc, ccp_node):
        if not isinstance(metacyc, MetaCyc):
            raise TypeError("The metacyc argument must be of the "
                            "MetaCyc class!")
        if not isinstance(ccp_node, Node):
            raise TypeError("The node argument must be of the "
                            "Node class or derived classes!")
        # TODO change ccp (?)
        edges_ccp = [e.source for e in metacyc.edges
                     if e.target == ccp_node and e.label == 'PART_OF']
        nodes = [n for n in metacyc.genes + metacyc.terminators +
                            metacyc.promoters + metacyc.BSs + metacyc.other_nodes
                 if hasattr(n, 'start') and hasattr(n, 'end') and n in edges_ccp]

        updated, several, notfound = [0]*3

        for node in nodes:
            location = self.identify_location(node.start, node.end)
            if location[2] is True:
                if len(location[0]):
                    node.start, node.end = location[0], location[1]
                    updated += 1
                else:
                    diff = [abs(l-node.start) for l in location[0]]
                    i = diff.index(min(diff))
                    node.start, node.end = location[0][i], location[1][i]
                    several += 1
                    # TODO LOGS
                    # TODO Label
            else:
                notfound += 1
                # TODO LOGS
                # TODO Label





