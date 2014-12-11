import warnings
import os
import re
import inspect as ins

###############################################################################


class CreateEdge():
    """
    An object of class CreateEdge.
    """
    def __init__(self, source, target, label=None):
        if not isinstance(source, Node):
            raise TypeError("The source argument must be of the Node class"
                            " or derived classes!")
        if not isinstance(target, Node):
            raise TypeError("The target argument must be of the Node class"
                            " or derived classes!")
        self.source = source
        self.target = target
        self.label = label

    def __repr__(self):
        return "An edge from %s to %s\n" \
               "Label: %s" %(self.source, self.target, self.label)

    def __str__(self):
        return "An edge from %s to %s\n" \
               "Label: %s" % (self.source, self.target, self.label)

    def __key(self):
        """
        The method returns a tuple of values stored in an edge, which is
        required for __eq__ and __hash__.
        """
        return tuple(self.__dict__.values())

    def __eq__(self, other):
        return self.__key() == other.__key()

    def __hash__(self):
        return hash(self.__key())

###############################################################################


class Node():
    """
    An object of class Node.
    """
    def __init__(self, source='MetaCyc'):
        self.labels = ':' + ':'.join(self.subclasses())
        self.comment = None
        self.source = source

    def __key(self):
        """
        The method returns a tuple of values stored in a node, which is
        required for __eq__ and __hash__.
        """
        return tuple(self.__dict__.values())

    def __eq__(self, other):
        return self.__key() == other.__key()

    def __hash__(self):
        return hash(self.__key())

    def __repr__(self):
        try:
            uid = self.uid
        except:
            uid = 'no uid specified'
        return "A node for %s class (uid: %s)" % (self.__class__.__name__, uid)

    def __str__(self):
        try:
            uid = self.uid
        except:
            uid = 'unknown'
        return "A node for %s class (uid: %s)" % (self.__class__.__name__, uid)

    def subclasses(self):
        """
        The method specifies a label attribute (string of class and all
        subclasses).
        """
        subclasses = [i.__name__ for i in ins.getmro(self.__class__) 
		      if i.__name__ != 'Node']
        return subclasses

###############################################################################


class BioEntity(Node):
    """
    Any molecular biology object with a name specified.
    """
    def __init__(self, name):
        if not isinstance(name, basestring):
            raise TypeError('The name argument must be a string!')
        Node.__init__(self)
        self.name = name

###############################################################################


class Feature(Node):
    """
    DNA regions with known biological features. Objects of class Feature have
    start, end and strand attributes.
    """
    def __init__(self, start, end, strand):
        if not isinstance(start, int):
            raise TypeError('The start argument must be an integer!')
        if not isinstance(end, int):
            raise TypeError('The end argument must be an integer!')
        if not isinstance(strand, basestring):
            raise TypeError('The strand argument must be a string!')
        if start < 0 or end < 0:
            warnings.warn("Check start and end positions! One of them "
                          "(or both) has a negative value!\n")
        if start > 0 and end > 0 and start > end:
            warnings.warn("Check start and end positions! The end "
                          "coordinate is smaller than the start "
                          "coordinate!\n")
        Node.__init__(self)
        self.start = start
        self.end = end
        self.strand = strand

###############################################################################

class Organism(Node):
    """
    An object of class Organism. It inherits all methods from the class Node.
    """
    def __init__(self, name):
        Node.__init__(self)
        self.name = name
        
###############################################################################

class Chromosome(BioEntity):
    """
    An object of class Chromosome. It inherits all methods from the class BioEntity.
    """
    def __init__(self, name, length, type):
        BioEntity.__init__(self, name)
        self.length = length
        self.type = type
        
###############################################################################

class Contig(BioEntity):
    """
    An object of class Contig. It inherits all methods from the class BioEntity.
    """
    def __init__(self, name, length, type):
        BioEntity.__init__(self, name)
        self.length = length
        self.type = type
        
###############################################################################

class Plasmid(BioEntity):
    """
    An object of class Plasmid. It inherits all methods from the class BioEntity.
    """
    def __init__(self, name, length, type):
        BioEntity.__init__(self, name)
        self.length = length
        self.type = type

###############################################################################

class Gene(BioEntity, Feature):
    """
    An object of class Gene. It inherits all methods from classes BioEntity
    and Feature.
    """
    def __init__(self, name, start, end, strand, uid=None, bcode=None,
                 product=None):
        BioEntity.__init__(self, name)
        Feature.__init__(self, start, end, strand)
        self.uid = uid
        self.bcode = bcode
        self.product = product

    def __repr__(self):
        return "%s gene" % self.name

    def __str__(self):
        return "Object of the class Gene:\n" \
               "name: %s\n " \
               "location: (%d, %d)\n" \
               "strand: %s" % (self.name, self.start, self.end, self.strand)

###############################################################################


class Term(Node):
    """
    An object of class Term. It inherits all methods from the Node class.
    """
    def __init__(self, text):
        if not isinstance(text, basestring):
            raise TypeError('The text argument must be a string!')
        Node.__init__(self)
        self.text = text
        self.source = None

###############################################################################


class XRef(Node):
    """
    An object of class XRef. It inherits all methods from the Node class.
    """
    def __init__(self, id):
        if not isinstance(id, basestring):
            raise TypeError('The id argument must be a string!')
        Node.__init__(self)
        self.id = id
        self.source = None

###############################################################################


class DB(Node):
    """
    An object of class DB. It inherits all methods from classes BioEntity and
    Feature.
    """
    def __init__(self, name, link=None):
        if not isinstance(name, basestring):
            raise TypeError('The name argument must be a string!')
        Node.__init__(self)
        self.name = name
        self.link = link
        self.source = None

###############################################################################


class RNA(BioEntity):
    """
    An object of class RNA. It inherits all methods from the class BioEntity.
    """
    def __init__(self, name, uid=None):
        BioEntity.__init__(self, name)
        self.uid = uid

###############################################################################


class rRNA(RNA):
    """
    An object of class rRNA. It inherits all methods from class RNA.
    """
    def __init__(self, name, uid=None):
        RNA.__init__(self, name)
        self.uid = uid

###############################################################################


class tRNA(RNA):
    """
    An object of class tRNA. It inherits all methods from class RNA.
    """
    def __init__(self, name, uid=None, type=None):
        RNA.__init__(self, name)
        self.uid = uid
        self.type = type

###############################################################################


class sRNA(RNA):
    """
    An object of class sRNA. It inherits all methods from class RNA.
    """
    def __init__(self, name, uid=None):
        RNA.__init__(self, name)
        self.uid = uid

###############################################################################


class Terminator(Feature):
    """
    An object of class Terminator. It inherits all methods from the class
    Feature.
    """
    def __init__(self, start, end, strand="unknown", uid=None):
        Feature.__init__(self, start, end, strand)
        self.uid = uid

###############################################################################


class Promoter(BioEntity, Feature):
    """
    An object of class Promoter. It inherits all methods from Feature and
    BioEntity classes.
    """
    def __init__(self, name, tss, start, end, strand="unknown", uid=None,
                 seq=None):
        BioEntity.__init__(self, name)
        Feature.__init__(self, start, end, strand)
        self.tss = tss
        self.uid = uid
        self.seq = seq

###############################################################################


class BS(Feature):
    """
    An object of class BS. It inherits all methods from the class Feature.
    """
    def __init__(self, start, end, strand="unknown", uid=None,
                 site_length=None):
        Feature.__init__(self, start, end, strand)
        self.uid = uid
        self.site_length = site_length

###############################################################################


class TU(BioEntity):
    """
    An object of class TU. It inherits all methods from the BioEntity class.
    """

    def __init__(self, name, uid=None):
        BioEntity.__init__(self, name)
        self.uid = uid

###############################################################################


class Peptide(BioEntity):
    """
    An object of class Peptide. It inherits all methods from the BioEntity
    class.
    """

    def __init__(self, name, uid=None, seq=None, molecular_weight_kd=None):
        BioEntity.__init__(self, name)
        self.uid = uid
        self.seq = seq
        self.molecular_weight_kd = molecular_weight_kd

###############################################################################


class Oligopeptide(Peptide):
    """
    An object of class Oligopeptide. It inherits all methods from the Peptide
    class.
    """
    def __init__(self, name, uid=None, seq=None, molecular_weight_kd=None):
        Peptide.__init__(self, name, uid, seq, molecular_weight_kd)

###############################################################################


class Polypeptide(Peptide):
    """
    An object of class Polypeptide. It inherits all methods from the Peptide
    class.
    """
    def __init__(self, name, uid=None, seq=None, molecular_weight_kd=None):
        Peptide.__init__(self, name, uid, seq, molecular_weight_kd)

###############################################################################


class Protein(BioEntity):
    """
    An object of class Protein. It inherits all methods from the BioEntity
    class.
    """
    def __init__(self, name):
        BioEntity.__init__(self, name)

###############################################################################


class Complex(BioEntity):
    """
    An object of class Complex. It inherits all methods from the BioEntity
    class.
    """
    def __init__(self, name, uid=None, molecular_weight_kd=None,
                 subunit_composition=None):
        BioEntity.__init__(self, name)
        self.uid = uid
        self.molecular_weight_kd = molecular_weight_kd
        self.subunit_composition = subunit_composition

###############################################################################


class Compound(Node):
    """
    An object of class Compound. It inherits all methods from the Node
    class.
    """

    def __init__(self, name, uid=None, chemical_formula=None, smiles=None,
                 type=None, molecular_weight=None):
        Node.__init__(self)
        self.uid = uid
        self.name = name
        self.chemical_formula = chemical_formula
        self.molecular_weight = molecular_weight
        self.smiles = smiles
        self.type = type

###############################################################################


class ProtFeature(Node):
    """
    An object of class ProtFeature. It inherits all methods from the Node
    class.
    """
    def __init__(self, uid=None, type=None, residues=None, source=None,
                 comment=None, alternate_sequence=None):
        Node.__init__(self)
        self.uid = uid
        self.type = type
        self.residues = residues
        self.source = source
        self.comment = comment
        self.alternate_sequence = alternate_sequence

###############################################################################


class SigmaFactor(Protein):
    """
    An object of class SigmaFactor. It inherits all methods from the Protein
    class.
    """
    def __init__(self, name):
        Protein.__init__(self, name)
        self.name = name

###############################################################################


class Enzyme(Protein):
    """
    An object of class Enzyme. It inherits all methods from the Protein class.
    """
    def __init__(self, name, uid=None):
        Protein.__init__(self, name)
        self.name = name
        self.uid = uid

###############################################################################


class Transporter(Protein):
    """
    An object of class Transporter. It inherits all methods from the Protein
    class.
    """
    def __init__(self, name, reaction=None):
        Protein.__init__(self, name)
        self.name = name
        self.reaction = reaction

###############################################################################


class RegulationEvent(Node):
    """
    An object of class RegulationEvent. It inherits all methods from the Node
    class.
    """
    def __init__(self, uid=None, comment=None):
        Node.__init__(self)
        self.uid = uid
        self.comment = comment

###############################################################################


class Binding(Node):
    """
    An object of class Binding. It inherits all methods from the Node class.
    """
    def __init__(self):
        Node.__init__(self)

###############################################################################


class Attenuation(RegulationEvent, Binding):
    """
    An object of class Attenuation. It inherits all methods from
    RegulationEvent and Binding classes.
    """
    def __init__(self, uid=None, comment=None, antiterminator_pos=None,
                 antiantiterminator_pos=None):
        RegulationEvent.__init__(self, uid, comment)
        Binding.__init__(self)
        self.uid = uid
        self.comment = comment
        self.antiterminator_pos = antiterminator_pos
        self.antiantiterminator_pos = antiantiterminator_pos

###############################################################################


class EnzymeRegulation(RegulationEvent):
    """
    An object of class EnzymeRegulation. It inherits all methods from the
    RegulationEvent.
    """
    def __init__(self, uid=None, comment=None, ki=None):
        RegulationEvent.__init__(self, uid, comment)
        self.uid = uid
        self.comment = comment
        self.ki = ki

###############################################################################


class TranscriptionRegulation(RegulationEvent, Binding):
    """
    An object of class TranscriptionRegulation. It inherits all methods from
    RegulationEvent and Binding classes.
    """
    def __init__(self, uid=None, comment=None):
        Binding.__init__(self)
        RegulationEvent.__init__(self, uid, comment)
        self.uid = uid
        self.comment = comment

###############################################################################


class TranslationRegulation(RegulationEvent, Binding):
    """
    An object of class TranslationRegulation. It inherits all methods from
    RegulationEvent and Binding classes.
    """
    def __init__(self, uid=None, comment=None):
        Binding.__init__(self)
        RegulationEvent.__init__(self, uid, comment)
        self.uid = uid
        self.comment = comment

###############################################################################


class Reaction(Node):
    """
    An object of class Reaction. It inherits all methods from the Node class.
    """
    def __init__(self, formula, uid=None, type=None):
        Node.__init__(self)
        self.formula = formula
        self.uid = uid
        self.type = type

###############################################################################


class Pathway(BioEntity):
    """
    An object of class Pathway. It inherits all methods from the BioEntity
    class.
    """
    def __init__(self, name, uid=None, type=None, reaction_layout=None):
        BioEntity.__init__(self, name)
        self.name = name
        self.uid = uid
        self.type = type
        self.reaction_layout = reaction_layout

###############################################################################


class Reactant(Node):
    """
    An object of class Reactant. It inherits all methods from the Node
    class.
    """
    def __init__(self, name, annotation, stoichiometric_coef=1):
        Node.__init__(self)
        self.name = name
        self.stoichiometric_coef = stoichiometric_coef
        self.annotation = annotation

###############################################################################


class Compartment(BioEntity):
    """
    An object of class Compartment. It inherits all methods from the BioEntity
    class.
    """
    def __init__(self, name, uid, organism):
        BioEntity.__init__(self, name)
        self.name = name
        self.uid = uid
        self.organism = organism

###############################################################################

class Unspecified(Node):
    """
    An object of class Unspecified. It inherits all methods from the Node
    class. It is used for all unclassified objects.
    """
    def __init__(self, name=None, uid=None):
        Node.__init__(self)
        self.name = name
        self.uid = uid

###############################################################################

#import doctest
#doctest.testfile("data_model_tests.txt")