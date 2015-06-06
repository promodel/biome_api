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
        new_start = self.new_ref.find(seq)
        if new_start == -1:
            return (False, False, False)
        else:
            return (new_start, new_start + len(seq) - 1, True)

    def update_metacyc(self, metacyc, ccp_node):
        if not isinstance(metacyc, MetaCyc):
            raise TypeError("The metacyc argument must be of the "
                            "MetaCyc class!")
        if not isinstance(ccp_node, Node):
            raise TypeError("The node argument must be of the "
                            "Node class or derived classes!")
        nodes = metacyc.genes + metacyc.terminators + metacyc.promoters +\
                metacyc.BSs + metacyc.other_nodes + metacyc.regulation_events
        edges_ccp = [e.source for e in metacyc.edges
                     if e.target == ccp_node and e.label == 'PART_OF']
        for node in nodes:
            if node in edges_ccp:
                location = [False, False, False]
                if hasattr(node, 'start') and hasattr(node, 'end'):
                    location = self.identify_location(node.start, node.end)
                if hasattr(node, 'antiterminator_start') and \
                        hasattr(node, 'antiterminator_end'):
                    location = self.identify_location(node.antiterminator_start,
                                                      node.antiterminator_end)
                if hasattr(node, 'antiantiterminator_start') and \
                        hasattr(node, 'antiantiterminator_end'):
                    location = self.identify_location(node.antiantiterminator_start,
                                                      node.antiantiterminator_end)
                if location[2] == True:
                    node.start, node.end = location[:2]




