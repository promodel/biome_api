import re
import os
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
        if not set(list(new_reference) + list(old_reference)) <= set(['A', 'T', 'G', 'C']):
            warnings.warn('Sequences have noncanonical nucleotides!')
        self.new_ref = new_reference
        self.old_ref = old_reference
        self.path = path

    def __repr__(self):
        return "_DatSet constructed for %s file" % self.filename

    def __str__(self):
        return "_DatSet constructed for %s file" % self.filename

