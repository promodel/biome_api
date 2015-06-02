from ...api import *
import warnings

class NewReference():
    """
    """
    def __init__(self, new_reference, old_reference, path='./'):
        new_reference = new_reference.upper()
        old_reference = old_reference.upper()
        if not isinstance(path, basestring):
            raise TypeError('The path argument must be a string!')
        if not isinstance(new_reference, basestring):
            raise TypeError('The new_reference argument must be a string!')
        if not isinstance(old_reference, basestring):
            raise TypeError('The old_reference argument must be a string!')
        if not isinstance(old_reference, basestring):
            raise TypeError('The old_reference argument must be a string!')
        if set(list(new_reference) + list(old_reference)) > set(['A', 'T', 'G', 'C']):
            warnings.warn('Sequences have noncanonical nucleotides!')
        self.new_ref = new_reference
        self.old_ref = old_reference
        self.path = path

    def __repr__(self):
        return "An object of NewReference class"

    def __str__(self):
        return "An object of NewReference class"

    def update_metacyc(self, metacyc, node):
        if not isinstance(metacyc, MetaCyc):
            raise TypeError("The metacyc argument must be of the "
                            "MetaCyc class!")
        if not isinstance(metacyc, Node):
            raise TypeError("The node argument must be of the "
                            "Node class or derived classes!")
        nodes = metacyc.genes + metacyc.terminators + metacyc.promoters +\
                metacyc.BSs + metacyc.other_nodes
        edges_ccp = [e.source for e in metacyc.edges
                     if e.target == node and e.label == 'PART_OF']
               



