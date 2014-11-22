from ...api import *
from Bio import SeqIO
from tabulate import tabulate
import networkx as nx

def items_list(objects):
    """
    The function removes empty elements from a list.
    """
    mylist = list(set(objects))
    if '' in mylist:
        mylist.remove('')
    return mylist

warnings.simplefilter('always', UserWarning)

def warning_on_one_line(message, category, filename, lineno, file=None,
                        line=None):
    """
    The function formats a warning message.
    """
    return ' %s:%s: %s:%s' % (filename, lineno, category.__name__, message)

warnings.formatwarning = warning_on_one_line

def make_name(name):
    """
    The function formats a string to be suitable name for storage.
    """
    return re.sub("[!,.'<>\[\];*-/+()\"]", "_", name).replace(' ', '_')

###############################################################################


class _DatSet():
    """
    A private class specially created to handle with .dat-files of MetaCyc.
    """
    def __init__(self, filename, path='./'):
        if not isinstance(path, basestring):
            raise TypeError('The path argument must be a string!')
        if not isinstance(filename, basestring):
            raise TypeError('The filename argument must be a string!')
        if not os.path.isdir(path):
            raise ValueError('The path does not exist!')
        self.filename = filename
        self.path = path
        self.data = None
        self.names = []

    def __repr__(self):
        return "_DatSet constructed for %s file" % self.filename

    def __str__(self):
        return "_DatSet constructed for %s file" % self.filename

    def readfile(self):
        """
        The methods works with dat-files format.
        """
        try:
            f = file('%s%s' % (self.path, self.filename), 'r')
            data = f.readlines()
            f.close()
            data = [line for line in data if line[0] != '#']
            chunks = {}
            chunk = _DatObject()
            i = 0
            for line in data:
                if line[:9] == 'UNIQUE-ID':
                    uid = line.replace("\n", '').split(" - ")[1]
                    self.names.append(uid)
                elif line[:2] == '//':
                    chunk.formatting()
                    chunks[uid] = chunk
                    chunk = _DatObject()
                elif line[:1] == '/':
                    pass
                else:
                    sp_line = line.replace('\n', '').split(" - ")
                    attr = make_name(sp_line[0])
                    if attr == "RIGHT" or attr == "LEFT":
                        sp_nextline = data[i+1].replace('\n', '').split(" - ")

                        # if it is the end of file
                        if len(data) != i + 2:
                            sp_nextnextline = data[i+2].replace('\n', '').split(" - ")
                        else:
                            sp_nextnextline = [None]

                        if sp_nextline[0] == '^COEFFICIENT':
                            sp_line[1] = '%s %s' % (sp_nextline[1], sp_line[1])

                        # information about compartment
                        if sp_nextline[0] == '^COMPARTMENT' :
                            cco = sp_nextline[1]
                            cco = make_name(cco)
                            if hasattr(chunk, cco):
                                attrval = str(getattr(chunk, cco))
                                attrval = '%s; %s' % (attrval, sp_line[1])
                                setattr(chunk, cco, attrval)
                            else:
                                setattr(chunk, cco, sp_line[1])
                        elif sp_nextnextline[0] == '^COMPARTMENT':
                            cco = sp_nextnextline[1]
                            cco = make_name(cco)
                            if hasattr(chunk, cco):
                                attrval = str(getattr(chunk, cco))
                                attrval = '%s; %s' % (attrval, sp_line[1])
                                setattr(chunk, cco, attrval)
                            else:
                                setattr(chunk, cco, sp_line[1])

                    # adding coefficients
                    if hasattr(chunk, attr):
                        attrval = str(getattr(chunk, attr))
                        attrval = '%s; %s' % (attrval, sp_line[1])
                        setattr(chunk, attr, attrval)
                    else:
                        setattr(chunk, attr, sp_line[1])

                i += 1
            self.data = chunks
            if len(chunks) == 0:
                warnings.warn("There are no data in the"
                              " %s file!\n" % self.filename)
        except:
            raise Exception()

###############################################################################


class _DatObject():
    """
    A private class for entries in .dat-files. Class methods test, format and
    upgrade entries and create a few specific links for nodes in a MetaCyc
    object based on the information in entries.
    """
    def __init__(self):
        """
        An entry could have different sets of attributes, that is why no
        attributes are predefined.
        """
        pass

    def dblink_format(self):
        """
        The method rewrite links to external Databases in more suitable
        for analysis format.
        """
        try:
            links = self.DBLINKS
            newlinks = []
            for link in links.split('; '):
                newlink = {}
                link = link.replace('(', '').replace(')', '').split(' ')
                if len(link) < 5:
                    continue
                newlink['DB'] = link[0]
                newlink['id'] = link[1].replace('"', '')
                newlink['DBnum'] = link[4]
                newlinks.append(newlink)
            self.DBLINKS = newlinks
        except:
            pass

    def location_format(self):
        """
        The method changes string type of location coordinates
        to integer type.
        """
        try:
            self.LEFT_END_POSITION = int(self.LEFT_END_POSITION)
            self.RIGHT_END_POSITION = int(self.RIGHT_END_POSITION)
        except:
            pass

    def tss_format(self):
        """
        The method changes string type of tss position to integer type.
        """
        try:
            self.ABSOLUTE_PLUS_1_POS = int(self.ABSOLUTE_PLUS_1_POS)
        except:
            pass

    def promoter_format(self):
        """
        The method tries to extract sequences from promoter entries
        comments and rewrite COMMENT attribute.
        """
        if hasattr(self, 'COMMENT'):
            try:
                seq = self.COMMENT.split("Promoter region "
                                         "sequence: ")[1].split(".")[0]
                seq = re.sub('[/{}]', "", seq)
                letters = [seq[i] for i in xrange(0, len(seq))]
                if set(letters) == set(['A', 'C', 'T', 'G']):
                    self.COMMENT = seq
                else:
                    delattr(self, "COMMENT")
            except:
                delattr(self, "COMMENT")
        else:
            pass

    def abs_center_format(self):
        """
        The method changes string type of binding site center position to
        integer/float type.
        """
        try:
            center = self.ABS_CENTER_POS.replace(" ", "")
            if center[-1] == '.':
                self.ABS_CENTER_POS = int(center[:-1])
            if center[-2:] == '.5':
                self.ABS_CENTER_POS = int(center[:-2]) + 0.5
            else:
                self.ABS_CENTER_POS = int(center)
        except:
            pass

    def mol_weight_format(self):
        """
        The method changes string type of molecular weight in KD to integer
        type.
        """
        try:
            mol_weight = self.MOLECULAR_WEIGHT_KD.replace(" ", "")
            self.MOLECULAR_WEIGHT_KD = float(mol_weight)
        except:
            pass

    def formatting(self):
        """
        The method collects all methods that are required for correct reading
        and formatting of data in .dat-file
        """
        self.dblink_format()
        self.location_format()
        self.mol_weight_format()
        if self.TYPES == 'Promoters':
            self.promoter_format()
            self.tss_format()
        if self.TYPES == 'DNA-Binding-Sites':
            self.abs_center_format()

    def links_to_db(self, node, metacyc):
        """
        The method takes as an input a Node subclass object and
        a MetaCyc database object and creates links to external databases
        for node. It forms (node)-[:EVIDENCE]->(xref)-[:LINK_TO]->(db)
        path and stores nodes and edges into the MetaCyc object.
        """
        if isinstance(node, Node):
            if isinstance(metacyc, MetaCyc):
                try:
                    for dblink in self.DBLINKS:
                        xref = XRef(dblink["id"])

                        # checking if there already exist a node with
                        # the same id
                        if xref in metacyc.xrefs:
                            i = metacyc.xrefs.index(xref)
                            xref = metacyc.xrefs[i]
                        else:
                            metacyc.xrefs.append(xref)

                        #checking if the database is already in self.dbs
                        db_obj = metacyc.db_checker(dblink["DB"])

                        # creating node -[EVIDENCE]-> xref -[LINK_TO]-> db_obj
                        # path
                        metacyc.edges.append(
                            CreateEdge(node, xref, 'EVIDENCE'))
                        metacyc.edges.append(
                            CreateEdge(xref, db_obj, 'LINK_TO'))
                except:
                    pass
            else:
                raise TypeError("The metacyc argument must "
                                "be of the MetaCyc class!")
        else:
            raise TypeError("The node argument must be of the Node "
                            "class or derived classes!")

    def links_to_synonyms(self, node, metacyc):
        """
        The method takes as an input a Node subclass object and
        a MetaCyc database object and creates links to node name synonyms.
        It forms (node)-[:HAS_NAME]->(term) link and stores nodes and
        edges into the MetaCyc object.
        """
        if isinstance(node, Node):
            if isinstance(metacyc, MetaCyc):
                try:
                    for synonym in self.SYNONYMS.split('; '):
                        synonym = synonym.replace(' ', '')
                        newterm = Term(synonym)

                        # checking if there already exist a node with
                        # the same text
                        if newterm in metacyc.terms:
                            i = metacyc.terms.index(newterm)
                            newterm = metacyc.terms[i]
                        else:
                            metacyc.terms.append(newterm)
                        metacyc.edges.append(
                            CreateEdge(node, newterm, 'HAS_NAME'))
                except:
                    pass
            else:
                raise TypeError("The metacyc argument must be of "
                                "the MetaCyc class!")
        else:
            raise TypeError("The node argument must be of the Node "
                            "class or derived classes!")

    def links_to_gene(self, node, metacyc):
        """
        The method takes as an input a Node subclass object and
        a MetaCyc database object and creates links from genes to node.
        It forms (gene)-[:ENCODES]->(node) link and stores nodes and
        edges into the MetaCyc object.
        """
        if isinstance(node, Node):
            if isinstance(metacyc, MetaCyc):
                try:
                    for s in self.GENE.split('; '):
                        gene = [g for g in metacyc.genes if (g.uid == s or
                                                             g.name == s)]
                        if len(gene) != 0:
                            metacyc.edges.append(
                                CreateEdge(gene[0], node, 'ENCODES'))
                except:
                    pass
            else:
                raise TypeError("The metacyc argument must be of "
                                "the MetaCyc class!")
        else:
            raise TypeError("The node argument must be of the Node "
                            "class or derived classes!")

    def links_to_components(self, node, metacyc):
        """
        The method takes as an input a Node subclass object and
        a MetaCyc database object and creates links from components to
        node. It forms (component)-[:PART_OF]->(node) link and stores
        nodes and edges into the MetaCyc object.
        """
        if isinstance(node, Node):
            if isinstance(metacyc, MetaCyc):
                try:
                    comps = self.COMPONENTS

                    # objects that could form complexes with others
                    db_objects = metacyc.rnas + metacyc.compounds\
                                 + metacyc.oligopeptides + metacyc.polypeptides
                    for comp in comps.split('; '):
                        component = [c for c in db_objects if c.uid == comp]

                        # checking if there is no component
                        if len(component) == 0:
                            print "The complex with %s is assigned before " \
                                  "its components (%s)!" % (node.uid, comp)
                            continue
                        metacyc.edges.append(
                            CreateEdge(component[0], node, 'PART_OF'))
                except:
                    pass
            else:
                raise TypeError("The metacyc argument must be of the MetaCyc "
                                "class!")
        else:
            raise TypeError("The node argument must be of the Node class or "
                            "derived classes!")

    def location_check(self, node):
        """
        The method checks if gene location is the same in genes.col
        and genes.dat files. After modification could be applied to
        any object of the class Feature.
        """
        if isinstance(node, Feature) and hasattr(node, 'uid'):
            if node.start == self.LEFT_END_POSITION and \
                            node.end == self.RIGHT_END_POSITION:
                pass
            else:
                warnings.warn("Location for a %s with unique id %s is "
                              "different in gene.dat and gene.col"
                              " files!" % (node.__class__.__name__, node.uid))
        else:
            raise TypeError("The node argument must be of the Feature class "
                            "or derived classes!")

    def name_check(self, node):
        """
        The method checks if gene name is the same in genes.col
        and genes.dat files. After modification could be applied to
        any object of the class BioEntity.
        """
        if isinstance(node, BioEntity) and hasattr(node, 'uid'):
            if (node.uid == node.name) or (node.name == self.COMMON_NAME):
                pass
            else:
                warnings.warn("Names for a %s with unique id %s are different "
                              "in gene.dat and gene.col "
                              "files!" % (node.__class__.__name__, node.uid))
        else:
            raise TypeError("The node argument must be of a BioEntity class or "
                            "a derived classes!")

    def attr_check(self, attrname, if_no_data=None):
        """
        The method helps to avoid errors in the creation of nodes when
        attribute is not specified in an entry. It checks if an attribute
        exists and if there is no data assign None or if_no_data value.
        """
        if isinstance(attrname, basestring):
            if hasattr(self, attrname):
                return getattr(self, attrname)
            else:
                return if_no_data
        else:
            raise TypeError("The attrname argument must be a string!")

    def feature_location(self):
        """
        The method returns formatted location of protein sequence features.
        """
        res = None
        if hasattr(self, "LEFT_END_POSITION"):
            res = '%d:%d' % (self.LEFT_END_POSITION, self.RIGHT_END_POSITION)
        elif hasattr(self, "RESIDUE_NUMBER"):
            res = self.RESIDUE_NUMBER
        return res

    def links_to_protein(self, node, metacyc):
        """
        The method takes as an input a Node subclass object and
        a MetaCyc database object and creates links from components to
        node. It forms (component)-[:PART_OF]->(node) link and stores
        nodes and edges into the MetaCyc object.
        """
        if isinstance(node, Node):
            if isinstance(metacyc, MetaCyc):
                try:
                    objects = metacyc.oligopeptides + metacyc.polypeptides\
                              + metacyc.complexes
                    objects += [p for p in metacyc.proteins
                                if p.__class__.__name__ not in ("SigmaFactor",
                                                                "Transporter")]
                    feature_of = self.FEATURE_OF.split('; ')
                    protein = [ob for ob in objects if ob.uid in feature_of]
                    if len(protein) != 0:
                        metacyc.edges.append(
                            CreateEdge(protein[0], node, 'HAS_FEATURE'))
                    else:
                        print self.FEATURE_OF
                except:
                    pass
            else:
                raise TypeError("The metacyc argument must be of the MetaCyc "
                                "class!")
        else:
            raise TypeError("The node argument must be of the Node class or "
                            "derived classes!")

    def links_to_regulated_entity(self, node, metacyc):
        """
        The method takes as an input a Node subclass object and a MetaCyc
        database object and creates links from node to regulated entity.
        It forms (node)-[:ACTIVATES/REPRESSES/UNKNOWN]->(regulated entity)
        link and stores nodes and edges into the MetaCyc object.
        """
        if isinstance(node, Node):
            if isinstance(metacyc, MetaCyc):
                try:
                    uid = self.REGULATED_ENTITY
                    enzymes = [p for p in metacyc.proteins
                               if p.__class__.__name__ == "Enzyme"]
                    objects = metacyc.promoters + enzymes + \
                              metacyc.terminators + metacyc.TUs
                    target = [t for t in objects if t.uid == uid]
                    if hasattr(self, "MODE"):
                        mode = self.MODE
                        if mode == '+':
                            edge_label = "ACTIVATES"
                        elif mode == '-':
                            edge_label = "REPRESSES"
                    else:
                        edge_label = 'UNKNOWN'
                    metacyc.edges.append(
                        CreateEdge(node, target[0], edge_label))
                except:
                    pass
            else:
                raise TypeError("The metacyc argument must be of the MetaCyc "
                                "class!")
        else:
            raise TypeError("The node argument must be of the Node class or "
                            "derived classes!")

    def links_to_regulator(self, node, metacyc):
        """
        The method takes as an input a Node subclass object and a MetaCyc
        database object and creates links from regulator to the node.
        It forms (regulator)-[:PARTICIPATES_IN]->(node) link and stores
        nodes and edges into the MetaCyc object.
        """
        if isinstance(node, Node):
            if isinstance(metacyc, MetaCyc):
                try:
                    uid = self.REGULATOR.replace("|", "")
                    objects = metacyc.compounds + metacyc.polypeptides + \
                              metacyc.oligopeptides + metacyc.rnas + \
                              metacyc.complexes + metacyc.other_nodes
                    source = [t for t in objects if t.uid == uid]
                    if len(source) == 0:
                        #print uid
                        # let's search in tRNAs
                        if uid[-5:] == "tRNAs":
                            rnas = [r for r in metacyc.rnas
                                    if r.__class__.__name__ == "tRNA" and r.type == uid]
                            for rna in rnas:
                                metacyc.edges.append(
                                    CreateEdge(rna, node, "PARTICIPATES_IN"))
                        # let's search in types of compounds
                        compounds = [c for c in metacyc.compounds
                                     if c.type == uid]
                        if len(compounds) != 0:
                            for compound in compounds:
                                metacyc.edges.append(
                                    CreateEdge(
                                        compound, node, "PARTICIPATES_IN"))
                    else:
                        metacyc.edges.append(
                            CreateEdge(source[0], node, "PARTICIPATES_IN"))
                except:
                    pass
            else:
                raise TypeError("The metacyc argument must be of the MetaCyc "
                                "class!")
        else:
            raise TypeError("The node argument must be of the Node class or "
                            "derived classes!")

    def links_to_enzymes(self, node, metacyc):
        """
        The method takes as an input a Node subclass object and a MetaCyc
        database object and creates links from an enzyme to the node. It forms
        (enzyme)-[:CATALYZES]->(node) link and stores nodes and edges into
        the MetaCyc object.
        """
        if isinstance(node, Node):
            #if isinstance(metacyc, MetaCyc):
            if metacyc.__class__.__name__ == 'MetaCyc':
                try:
                    enzymes = [p for p in metacyc.proteins
                               if p.__class__.__name__ == "Enzyme"]
                    for item in self.ENZYMATIC_REACTION.split('; '):
                        enzyme = [e for e in enzymes if e.uid == item][0]
                        metacyc.edges.append(
                            CreateEdge(enzyme, node, 'CATALYZES'))
                except:
                    pass
            else:
                raise TypeError("The metacyc argument must be of the MetaCyc "
                                "class!")
        else:
            raise TypeError("The node argument must be of the Node class or "
                            "derived classes!")

    def links_to_reactions(self, node, metacyc):
        """
        The method takes as an input a Node subclass object and a MetaCyc
        database object and creates links from a reaction to the node. It forms
        (reaction/pathway)-[:PART_OF]->(node) link and stores nodes and edges
        into the MetaCyc object.
        """
        if isinstance(node, Node):
            if isinstance(metacyc, MetaCyc):
                try:
                    reacts_paths = metacyc.reactions + metacyc.pathways
                    reactions = [r for r in reacts_paths
                                 if r.uid in self.REACTION_LIST.split('; ')]
                    for reaction in reactions:
                        metacyc.edges.append(
                            CreateEdge(reaction, node, 'PART_OF'))
                except:
                    pass
            else:
                raise TypeError("The metacyc argument must be of the MetaCyc "
                                "class!")
        else:
            raise TypeError("The node argument must be of the Node class or "
                            "derived classes!")

    def match_compartments(self):
        """
        """
        cco = {
            'CCO-EXTRACELLULAR-CCO-IN': ('Extracellular space', 'Cytosol'),
            'CCO-PM-BAC-POS': ('Cell wall', 'Cytosol'),
            'CCO-EXTRACELLULAR-CCO-CYTOSOL': ('Extracellular space',
                                              'Cytosol'),
            'CCO-CE-BAC-POS': ('Cell envelope',),
            'CCO-CYTOSOL': ('Cytosol',),
            'CCO-PM-BAC-NEG': ('Periplasmic space', 'Cytosol')}

        reactants = self.LEFT.split('; ') + self.RIGHT.split('; ')

        # if no data about compartments, the default compartment is
        # cytosol
        if not hasattr(self, 'RXN_LOCATIONS') and \
                not hasattr(self, '^COMPARTMENT'):
            return [('CCO-CYTOSOL', 'Cytosol', reactant)
                    for reactant in reactants]

        # if reaction substrates and products are located in one
        # compartment and the compartment is specified in the
        # RXN_LOCATIONS slot
        elif hasattr(self, 'RXN_LOCATIONS') and \
                not hasattr(self, '^COMPARTMENT'):
            if getattr(self, 'RXN_LOCATIONS') in cco.keys():
                compartment = getattr(self, 'RXN_LOCATIONS')
                return [(compartment, cco[compartment][0], reactant)
                        for reactant in reactants]
            else:
                return [('UNKNOWN', 'Unknown', reactant)
                        for reactant in reactants]

        # if reaction substrates and products are located in different
        # compartments and the compartments are specified in the
        # ^COMPARTMENT slot
        elif not hasattr(self, 'RXN_LOCATIONS') and \
                hasattr(self, '^COMPARTMENT'):
            res = []
            compartments = list(set(getattr(self, '^COMPARTMENT').split('; ')))
            for compartment in compartments:
                reactants = getattr(self, make_name(compartment)).split('; ')
                if compartment in cco.keys():
                    cco_loc = cco[compartment][0]
                    res1 = [(compartment, cco_loc, reactant)
                            for reactant in reactants]
                else:
                    res1 = [('UNKNOWN', 'Unknown', reactant)
                            for reactant in reactants]
                res.extend(res1)
            return res

        # if reaction substrates and products are located in different
        # compartments and the compartments are specified in the
        # RXN_LOCATIONS slot
        elif hasattr(self, 'RXN_LOCATIONS') and hasattr(self, '^COMPARTMENT')\
                and (hasattr(self, 'CCO_IN') or hasattr(self, 'CCO_OUT')):
            compartment = getattr(self, 'RXN_LOCATIONS')
            try:
                cco_in = cco[compartment][0]
                res_in = [(compartment, cco_in, reactant)
                          for reactant in self.CCO_IN.split('; ')]
            except:
                res_in = []
            try:
                cco_out = cco[compartment][1]
                res_out = [(compartment, cco_out, reactant)
                           for reactant in self.CCO_OUT.split('; ')]
            except:
                res_out = []

            cco_in_out = []
            if hasattr(self, 'CCO_IN'):
                cco_in_out.append(self.CCO_IN)
            if hasattr(self, 'CCO_OUT'):
                cco_in_out.append(self.CCO_OUT)

            res_unknowns = [('UNKNOWN', 'Unknown', reactant)
                            for reactant in reactants
                            if reactant not in cco_in_out]

            res = res_in + res_out + res_unknowns
            return res

        else:
            raise Exception('Unexpected error in match_compartments!')

    def make_formula(self):
        left = self.LEFT.split('; ')
        right = self.RIGHT.split('; ')

        if hasattr(self, "REACTION_DIRECTION"):
            if self.REACTION_DIRECTION == "LEFT-TO-RIGHT" or \
                            self.REACTION_DIRECTION == "PHYSIOL-LEFT-TO-RIGHT":
                direction = "-->"
            elif self.REACTION_DIRECTION == "RIGHT-TO-LEFT" or \
                            self.REACTION_DIRECTION == "PHYSIOL-RIGHT-TO-LEFT":
                direction = "<--"
            elif self.REACTION_DIRECTION == "REVERSIBLE":
                direction = "<-->"
            else:
                direction = "--"
        else:
            direction = "--"
        formula = "%s %s %s" % (" + ".join(left), direction,
                                " + ".join(right))
        return formula

    def links_to_reactants(self, node, metacyc):
        """
        The method is looking for reaction reactants in MetaCyc object nodes.
        If there are no reactants the method creates them as Unspecified
        objects. As output it returns a tuple with a list of reactants and
        a number: 0 (if reactants were in MetaCyc object nodes) or a number
        of nodes it has created.
        """
        if isinstance(node, Node):
            if isinstance(metacyc, MetaCyc):
                objects = metacyc.compounds + metacyc.complexes + \
                          metacyc.polypeptides + metacyc.oligopeptides + \
                          metacyc.other_nodes
                groups = [rna for rna in metacyc.rnas if
                          rna.__class__.__name__ == 'tRNA'] + metacyc.compounds
                comps = self.match_compartments()

                for comp in comps:
                    # separating reactant name and numerical stoichiometric
                    # coefficient
                    reagent = comp[2]
                    item = re.sub(r'\d\s', '', reagent)

                    try:
                        num = re.match('\d+', reagent).group(0)
                    except:
                        num = 1

                    # if stoichiometric coefficient is equal to n
                        if reagent[:2] == 'n ':
                            substance = reagent[2:]
                            num = 'n'

                    reagents = [i for i in objects
                                       if i.uid == reagent or i.name == reagent]


                    # searching in names of object groups and classes
                    if len(reagents) == 0:
                        reagents = [i for i in groups
                                    if reagent in i.type.split('; ')]

                    # creating nodes for those reactant that were not
                    # found
                    if len(reagents) == 0:
                        notfound = Unspecified(name=item.replace('|', ''))
                        if notfound in metacyc.other_nodes:
                            i = metacyc.other_nodes.index(notfound)
                            reagents = [metacyc.other_nodes[i]]
                        else:
                            metacyc.other_nodes.append(notfound)
                            reagents = [notfound]

                    for reagent in reagents:
                        compartment = [item for item in metacyc.compartments
                                       if item.name == comp[1]]
                        if len(compartment) == 0:
                            compartment = [item for item in metacyc.compartments
                                       if item.name == 'Unknown']
                        reactant_name = '%s [%s]' % (reagent.name, comp[1])
                        annotation = '%s_%s_%s' % (reagent.name, comp[1],
                                                   metacyc.organism)
                        reactant = Reactant(stoichiometric_coef=num,
                                            name=reactant_name,
                                            annotation = annotation)
                        metacyc.reactants.append(reactant)
                        metacyc.edges.append(
                            CreateEdge(reagent, reactant, 'IS_A'))
                        metacyc.edges.append(
                            CreateEdge(reactant, node, 'PARTICIPATES_IN'))
                        metacyc.edges.append(
                            CreateEdge(reactant, compartment[0], 'LOCATES_IN'))
            else:
                raise TypeError("The metacyc argument must be of the MetaCyc "
                                    "class!")
        else:
            raise TypeError("The node argument must be of the Node class or "
                            "derived classes!")

    def reactions_order(self, metacyc):
        """
        The method tries to create chains of molecular reaction. It forms
        (reaction1) -[:NEXT]-> (reaction2) -[:NEXT]-> (reaction3) paths.
        """
        if isinstance(metacyc, MetaCyc):
            try:
                for pred in self.PREDECESSORS.split('; '):
                    pred = re.sub("()", "", pred).split(" ")
                    reaction1 = [r for r in metacyc.reactions
                                 if r.uid == pred[1]][0]
                    reaction2 = [r for r in metacyc.reactions
                                 if r.uid == pred[0]][0]
                    check_edge = [e for e in metacyc.edges
                                  if e.source == reaction1 and e.target == reaction2]
                    if len(check_edge) == 0:
                        metacyc.edges.append(
                            CreateEdge(reaction1, reaction2, 'NEXT'))
                    else:
                        pass
            except:
                pass
        else:
            raise TypeError("The metacyc argument must be of the MetaCyc "
                            "class!")

###############################################################################


class MetaCyc():
    """
    Class for data taken from the MetaCyc database.
    """
    def __init__(self, path=''):
        if not isinstance(path, basestring):
            raise TypeError('The path argument must be a string!')
        if not os.path.isdir(path):
            raise ValueError('The path does not exist!')
        self.path = path
        self.organism_name = None
        self.orgid = None
        self.version = None
        self.release = None
        self.organism = []
        self.ccp = []
        self.genes = []
        self.edges = []
        self.xrefs = []
        self.dbs = []
        self.terms = []
        self.rnas = []
        self.terminators = []
        self.promoters = []
        self.BSs = []
        self.TUs = []
        self.compounds = []
        self.polypeptides = []
        self.oligopeptides = []
        self.proteins = []
        self.complexes = []
        self.protfeatures = []
        self.regulation_events = []
        self.reactions = []
        self.reactants = []
        self.compartments = []
        self.pathways = []
        self.other_nodes = []

    def __repr__(self):
        if self.organism != None and self.version != None\
                and self.release != None:
            return "MetaCyc object for %s: version %s, data of " \
                   "release %s" % (self.organism, self.version, self.release)
        else:
            return "MetaCyc object"

    def __str__(self):
        if self.organism != None and self.version != None \
                and self.release != None:
            return "MetaCyc object for %s: version %s, data of " \
                   "release %s" % (self.organism, self.version, self.release)
        else:
            return "MetaCyc object"

    def _read_dat(self, filename):
        """
        The method tries to read .dat-file
        """
        try:
            # reading the file
            datfile = _DatSet(filename, self.path)
            datfile.readfile()
            return datfile
        except:
            print "There is no %s file in the database or it has wrong " \
                  "format! Let's skip it..." % filename
            return None

    def add_edge(self, source, target, label):
        """
        The method creates an edge only if it is not presented in the
        MetaCyc object edges attribute.
        """
        edge = CreateEdge(source, target, label)
        if edge not in self.edges:
            self.edges.append(edge)

    def _location(self, start, end, trans_dir=None):
        """
        Tne method takes start and end positions as an input and returns a
        correct list of the start, the end, the strand for an object based
        on input values.
        """
        if not isinstance(start, int):
            raise TypeError('The start argument must be an integer!')
        if not isinstance(end, int):
            raise TypeError('The end argument must be an integer!')
	if trans_dir not in ['+', '-', None]:
            raise ValueError('The trans_dir argument must be an integer!')
	if trans_dir == '+':
	  return [start, end, 'forward']
	elif trans_dir == '-':
	  return [start, end, 'reverse']
	else:
	  if start < end:
	      return [start, end, 'forward']
	  elif start > end:
	      return [end, start, 'reverse']
	  else:
	      return [start, end, 'unknown']
	

    def _set_version(self):
        """
        A private method that assigns DB version to an object of the class
        MetaCyc.
        """
        try:
            f = file(self.path + "version.dat", 'r')
            data = f.readlines()
            f.close()
            for line in data:
                if line[0] != ';':
                    chunks = line.replace('\n', '').split('\t')
                    if chunks[0] == 'ORGID':
                        self.orgid = chunks[1]
                    elif chunks[0] == 'ORGANISM':
                        self.organism_name = chunks[1]
                    elif chunks[0] == 'VERSION':
                        self.version = chunks[1]
                    elif chunks[0] == 'RELEASE-DATE':
                        self.release = chunks[1]
            print "Information about version and release has been set!"
        except:
            print 'There is no information about the database version!'
            pass

    def extract_data(self):
        """
        Tne method uses a number of methods for data extraction.
        """
        # Setting version information
        self._set_version()

        # Everything about genes; methods create nodes for genes, terms, xrefs,
        # dbs and edges between them
        self.genes_col()
        self.genes_dat()
        self.gene_links()

        # RNAs
        self.rnas_dat()

        # Terminators
        self.terminators_dat()

        # Promoters
        self.promoters_dat()

        # BSs
        self.dnabindsites_dat()

        # Transcription units
        self.transunits_dat()

        # Compounds
        self.compounds_dat()

        # Proteins, peptides and their complexes
        self.proteins_dat()
        self.protseq_fsa()
        self.protcplxs_col()
        self.protein_features_dat()
        self.sigma_factors()
        self.enzrxns_dat()
        self.transporters_col()

        # Regulation
        self.regulation_dat()

        # Pathways
        self.create_compartments()
        self.reactions_dat()
        self.pathways_dat()
        self.pathways_col()

    def db_checker(self, name):
        """
        The method checks if there is already a node for a particular database.
        """
        if not isinstance(name, basestring):
            raise TypeError('The name argument must be a string!')
        if len(self.dbs) == 0:
            db_obj = DB(name)
            self.dbs.append(db_obj)
        else:
            db_obj = [d for d in self.dbs if d.name == name]
            if len(db_obj) == 1:
                db_obj = db_obj[0]
            elif len(db_obj) == 0:
                db_obj = DB(name)
                self.dbs.append(db_obj)
            else:
                print "Unexpected error!"
        return db_obj

    def name_to_terms(self, node):
        if isinstance(node, Node):
            if hasattr(node, 'name'):
                term = Term(text = node.name)
                self.terms.append(term)
                self.edges.append(CreateEdge(node, term, 'HAS_NAME'))
            else:
                pass
        else:
            raise TypeError("The node argument must be of the Node class or "
                            "derived classes!")

    def genes_col(self):
        """
        Data extraction from genes.col file. Only writing (no upgrading)
        is possible!
        """
        try:
            f = file(self.path + "genes.col", 'r')
            data = f.readlines()
            f.close()

            # WRITING de novo
            if len(self.genes) == 0:
                for line in data:
                    if line[0] != '#' and line[:9] != 'UNIQUE-ID':
                        chunks = line.replace('\n', '').split('\t')

                        # creating formatted coordinates for a gene
                        location = self._location(
                            int(chunks[5]), int(chunks[6]))

                        gene = Gene(uid=chunks[0], name=chunks[1],
                                    start=location[0], end=location[1],
                                    strand=location[2], product=chunks[2])
                        self.genes.append(gene)
                        self.name_to_terms(gene)
                print "A list with %d genes has been " \
                      "created!" % len(self.genes)

            # UPGRADING
            else:
                pass
        except:
            print "There is no genes.col file in the database or it has " \
                  "wrong format! Let's skip it..."

    def genes_dat(self):
        """
        The method stores data about genes from genes.dat file to MetaCyc
        object or upgrades already existing information about genes.
        """
        datfile = self._read_dat('genes.dat')
        if datfile is not None:
            # WRITING de novo
            if len(self.genes) == 0 and datfile.data is not None:

                # creating Gene objects; not all attributes are presented for
                # each _DatObject so it's important to use attr_check method
                for uid in datfile.names:
                    obj = datfile.data[uid]
                    location = self._location(
                            obj.attr_check("LEFT_END_POSITION"), 
                            obj.attr_check("RIGHT_END_POSITION"),
                            obj.attr_check("TRANSCRIPTION-DIRECTION"))
                    gene = Gene(uid=uid,
                                name=obj.attr_check("COMMON_NAME", uid),
                                start=location[0],
                                end=location[1],
                                strand=location[2],
                                product=obj.attr_check("PRODUCT"))
                    self.genes.append(gene)
                    self.name_to_terms(gene)

                    # creating Terms for gene name synonyms
                    obj.links_to_synonyms(gene, self)

                    # creating XRef and DB nodes for gene name dblinks
                    obj.links_to_db(gene, self)

                print "A list with %d genes have been " \
                      "created!" % len(self.genes)

            # UPGRADING of self.genes
            elif len(self.genes) != 0 and datfile.data is not None:
                for gene in self.genes:
                    uid = gene.uid
                    try:
                        obj = datfile.data[uid]

                        # let's check that genes with the same id have
                        # the same name
                        #obj.name_check(gene)

                        # let's check that genes with the same id have
                        # the same location
                        obj.location_check(gene)

                        # creating Terms for gene name synonyms
                        obj.links_to_synonyms(gene, self)

                        # creating XRef and DB nodes for gene name dblinks
                        obj.links_to_db(gene, self)
                    except:
                        pass
            else:
                raise StandardError("Something wrong has happened!")

    def gene_links(self):
        """
        The method creates links from genes to Uniprot entries,
        It forms (gene) -[:EVIDENCE]-> (xref) -[:LINK_TO]-> (DB:Uniprot) path.
        """
        if len(self.genes) != 0:
            try:
                f = file(self.path + "gene-links.dat", 'r')
                data = f.readlines()
                f.close()
                for line in data:
                    if line[0] == "#":
                        continue
                    chunks = line.replace('\n', '').split('\t')
                    if chunks[1] == '':
                        continue
                    gene = [g for g in self.genes if g.uid == chunks[0]]
                    if len(gene) == 0:
                        continue

                    xref = XRef(chunks[1])
                    self.xrefs.append(xref)
                    self.edges.append(CreateEdge(gene[0], xref, 'EVIDENCE'))
                    db_obj = self.db_checker("UniProt")
                    self.edges.append(CreateEdge(xref, db_obj, 'LINK_TO'))
            except:
                print "There is no gene-links.dat file in the database or " \
                      "it has wrong format! Let's skip it..."
        else:
            raise TypeError("There is no information about genes in the "
                            "database! Call genes_col or genes_dat method "
                            "first!")

    def rnas_dat(self):
        """
        The method creates three types of RNA nodes (rRNA, tRNA, sRNA) and
        links to their modified forms.
        """
        datfile = self._read_dat('rnas.dat')
        if datfile is not None:
            notfound = 0
            for uid in datfile.names:
                obj = datfile.data[uid]

                # we skip all modified rnas records
                if hasattr(obj, "UNMODIFIED_FORM"):
                    continue
                if obj.TYPES[:6] == 'Charge':
                    continue

                # picking gene that encodes the rna
                gene = [g for g in self.genes
                        if g.uid == obj.GENE or g.name == obj.GENE]

                # creating rna object (3 types)
                if len(gene) == 1:
                    gene = gene[0]
                    if obj.TYPES[-5:] == '-RNAs':
                        rna = sRNA(uid=uid,
                                   name=obj.attr_check("COMMON_NAME", uid))
                    if obj.TYPES[-5:] == 'tRNAs':
                        rna = tRNA(uid=uid,
                                   name=obj.attr_check("COMMON_NAME", uid),
                                   type=obj.attr_check("TYPES"))
                    if obj.TYPES[-5:] == 'rRNAs':
                        rna = rRNA(uid=uid,
                                   name=obj.attr_check("COMMON_NAME", uid))
                else:
                    notfound += 1
                    continue

                self.rnas.append(rna)
                self.name_to_terms(rna)

                # creating Terms for RNAs name synonyms
                obj.links_to_synonyms(rna, self)

                # creating edges gene -[ENCODES]-> rna
                self.edges.append(CreateEdge(gene, rna, 'ENCODES'))

                # creating a node for rna modification (if it exists)
                if hasattr(obj, "MODIFIED_FORM"):
                    uid = obj.MODIFIED_FORM.replace("|", "")
                    mod_obj = datfile.data[uid]
                    modification = Unspecified(
                        name=mod_obj.attr_check("COMMON_NAME", uid), uid=uid)
                    self.other_nodes.append(modification)
                    self.name_to_terms(modification)

                    # creating edges rna -[FORMS]-> modification
                    self.edges.append(CreateEdge(rna, modification, 'FORMS'))

            print "A list with %d RNAs has been created!" % len(self.rnas)
            print "No genes in the database for %d RNAs." % notfound

    def terminators_dat(self):
        """
        The method creates nodes for terminators.
        """
        datfile = self._read_dat('terminators.dat')
        if datfile is not None:
            unmapped =0
            for uid in datfile.names:
                obj = datfile.data[uid]
                if hasattr(obj, "UNMAPPED_COMPONENT_OF"):
                    unmapped += 1
                    continue
                ter = Terminator(uid=uid, start=obj.LEFT_END_POSITION,
                                 end=obj.RIGHT_END_POSITION)
                self.terminators.append(ter)
            print "A list with %d terminators has been created!\n" \
                  "There were %d unmapped elements, they were " \
                  "skipped..." % (len(self.terminators), unmapped)

    def promoters_dat(self):
        """
        The method creates nodes for promoters.
        """
        datfile = self._read_dat('promoters.dat')
        if datfile is not None:
            unmapped = 0
            for uid in datfile.names:
                obj = datfile.data[uid]
                if not hasattr(obj, "ABSOLUTE_PLUS_1_POS"):
                    unmapped += 1
                    continue
                pro = Promoter(uid=uid,
                               name=obj.attr_check("COMMON_NAME", uid),
                               start=obj.ABSOLUTE_PLUS_1_POS,
                               end=obj.ABSOLUTE_PLUS_1_POS,
                               tss=obj.ABSOLUTE_PLUS_1_POS,
                               strand="unknown", seq=None)
                self.promoters.append(pro)
                self.name_to_terms(pro)

            print "A list with %d promoters has been created!\n" \
                  "There were %d unmapped elements, they were " \
                  "skipped..." % (len(self.promoters), unmapped)

    def dnabindsites_dat(self):
        """
        The method creates nodes for binding sites.
        """
        datfile = self._read_dat('dnabindsites.dat')
        if datfile is not None:
            unmapped = 0
            for uid in datfile.names:
                obj = datfile.data[uid]

                # checking if binding site location is specified in an entry
                if not hasattr(obj, "ABS_CENTER_POS"):
                    unmapped += 1
                    continue

                # formatting binding site location
                if isinstance(obj.ABS_CENTER_POS, float):
                    bsite = BS(uid=uid, start=int(obj.ABS_CENTER_POS - 0.5),
                               end=int(obj.ABS_CENTER_POS + 0.5),
                               site_length=obj.attr_check("SITE_LENGTH"))
                    self.BSs.append(bsite)
                else:
                    bsite = BS(uid=uid, start=obj.ABS_CENTER_POS,
                               end=obj.ABS_CENTER_POS,
                               site_length=obj.attr_check("SITE_LENGTH"))
                    self.BSs.append(bsite)
            print "A list with %d BSs has been created!\n" \
                  "There were %d unmapped BSs, they were " \
                  "skipped..." % (len(self.BSs), unmapped)

    def create_organism(self):
        """
        The method creates an Organism node
        """
        self.organism.append(Organism(name=self.organism_name))

    def create_ccp(self):
        """
        The method reads data in .nt-file and creates chromosomes, contigs
        or plasmids.
        """
        try:
        # reading genome sequence from .nt-file
            for mcfile in os.listdir(self.path):
                if mcfile.endswith(".nt"):
                    filename = mcfile
            f = open(self.path + filename, "rU")
            records = list(SeqIO.parse(f, "fasta"))
            f.close()
        except:
            print "There is no .nt-file!"
            return None

        # Creating chromosomes, contigs or plasmids
        for record in records:
            name = record.description.split('|')[-1]
            length = len(record.seq)
            record_id = record.name.split('|')[-1]
            print record_id, length, name
            if ('complete genome' in record.description or \
                        'complete sequence' in record.description) and \
                            'lasmid' not in record.description:
                ccp_obj = Chromosome(name=name, length=length,
                                     accesion=record_id, type='unknown')
            elif 'ontig' in record.description:
                ccp_obj = Contig(name=name, length=length,
                                 accesion=record_id, type='unknown')
            elif 'lasmid' in record.description:
                ccp_obj = Plasmid(name=name, length=length,
                                  accesion=record_id, type='unknown')
            else:
                raise UserWarning('Unknown genome element')

            self.ccp.append(ccp_obj)
            self.edges.append(
                CreateEdge(ccp_obj, self.organism, 'PART_OF'))

    def transunits_dat(self):
        """
        The method creates nodes for transcription units and create links to
        elements they contain. It also rewrite strand attributes in the
        elements and specify promoter sequence, makes a few tests for elements
        location.
        """
        datfile = self._read_dat('transunits.dat')
        if datfile is not None:
            notcomplete = 0
            nocomps = 0
            elements = self.promoters + self.BSs + self.genes\
                       + self.terminators
            for uid in datfile.names:
                obj = datfile.data[uid]
                tu_obj = TU(name=obj.attr_check("COMMON_NAME", uid), uid=uid)

                comps = [c for c in elements \
                         if c.uid in obj.COMPONENTS.split("; ")]
                if len(comps) == 0:
                    nocomps += 1
                    continue
                if len(comps) != len(obj.COMPONENTS.split("; ")):
                    notcomplete += 1

                self.TUs.append(tu_obj)
                self.name_to_terms(tu_obj)

                genes = [c for c in comps if c.__class__.__name__ == "Gene"]

                # if there are no genes
                if len(genes) == 0:
                    for comp in comps:
                        self.edges.append(
                            CreateEdge(tu_obj, comp, 'CONTAINS'))
                    continue

                # checking genes location and strand values
                test = _GenesTest(genes)

                # creating edges and rewriting strand values to promoters,
                # terminators and BSs
                for comp in comps:
                    comp_class = comp.__class__.__name__
                    if comp_class == "Gene":
                        self.edges.append(
                            CreateEdge(tu_obj, comp, 'CONTAINS'))
                    else:
                        comp.strand = test.strand
                        # rewriting coordinates and sequences for promoters
                        if comp_class == "Promoter":
                            if comp.strand == 'forward':
                                comp.start = comp.start - 60
                                comp.end = comp.end + 15
                            elif comp.strand == 'reverse':
                                comp.start = comp.start - 15
                                comp.end = comp.end + 60
                            if "genome" in locals():
                                if comp.start < 0:
                                    comp.start = 0
                                elif comp.start > len(genome):
                                    comp.start = len(genome)
                                if comp.end < 0:
                                    comp.end = 0
                                elif comp.end > len(genome):
                                    comp.end = len(genome)
                                seq = genome[comp.start + 1 : comp.end + 1]
                                comp.seq = seq
                        self.edges.append(
                            CreateEdge(tu_obj, comp, 'CONTAINS'))
            print "A list with %d transcription units has been created!\n" \
                  "There were %d incomplete transcription units...\n" \
                  "No information about %d transcription " \
                  "units." % (len(self.TUs), notcomplete, nocomps)

    def compounds_dat(self):
        """
        The method creates nodes for compounds. It also makes links to
        external Databases and compound name synonyms.
        """
        datfile = self._read_dat('compounds.dat')
        if datfile is not None:
            for uid in datfile.names:
                obj = datfile.data[uid]
                compound = Compound(uid=uid,
                                    name=obj.attr_check("COMMON_NAME", uid),
                                    chemical_formula=obj.attr_check("CHEMICAL_FORMULA"),
                                    molecular_weight=obj.attr_check("MOLECULAR_WEIGHT"),
                                    smiles=obj.attr_check("SMILES"),
                                    type=obj.attr_check("TYPES"))
                self.compounds.append(compound)
                self.name_to_terms(compound)

                # creating Terms for compounds name synonyms
                obj.links_to_synonyms(compound, self)

                # creating XRef and DB nodes for compounds dblinks
                obj.links_to_db(compound, self)

            print "A list with %d compounds has been " \
                  "created!" % len(self.compounds)

    def proteins_dat(self):
        """
        The method creates nodes for various peptides and proteins. It also
        makes links to external Databases, peptide/protein name synonyms,
        genes and their modified forms.
        """
        datfile = self._read_dat('proteins.dat')
        if datfile is not None:
            # There are 2 types of peptides: oligopeptides, polypeptides
            oligo = ['OLIGOPEPTIDES', 'Cycic-Dipeptides',
                     'Dipeptides-Of-D-Amino-Acids', 'DIPEPTIDES',
                     'Tetrapeptides', 'TRIPEPTIDES']
            poly = ['Polypeptides', 'apo-ACP', 'Reduced-2Fe-2S-Ferredoxins',
                    'SpoIVB-peptidase-precursors', 'All-ACPs']
            complexes = []
            for uid in datfile.names:
                obj = datfile.data[uid]

                # we skip all modified forms, as unmodified peptides entries
                # have information about them
                if hasattr(obj, 'UNMODIFIED_FORM'):
                    continue

                # we also skip complexes, just to be sure that later all
                # elements for them will exist
                if uid[:4] == 'CPLX':
                    complexes.append(uid)
                    continue

                types = obj.TYPES.split('; ')

                # Polypeptides
                if len(set(types + poly)) != len(poly) + len(types):
                    peptide = Polypeptide(uid=uid,
                                          name=obj.attr_check("COMMON_NAME", uid),
                                          molecular_weight_kd=obj.attr_check("MOLECULAR_WEIGHT_KD"))
                    self.polypeptides.append(peptide)

                # Oligopeptides
                elif len(set(types + oligo)) != len(oligo) + len(types):
                    peptide = Oligopeptide(uid=uid,
                                           name=obj.attr_check("COMMON_NAME", uid),
                                           molecular_weight_kd=obj.attr_check("MOLECULAR_WEIGHT_KD"))
                    self.oligopeptides.append(peptide)
                else:
                    warnings.warn("Unexpected peptide types! "
                                  "Let's skip it... \nThe object uid %s" % uid)
                    continue
                self.name_to_terms(peptide)
                
                # creating Terms for peptide name synonyms
                obj.links_to_synonyms(peptide, self)

                # creating XRef and DB nodes for peptide dblinks
                obj.links_to_db(peptide, self)

                # creating XRef and DB node for GO-terms
                if hasattr(obj, 'GO_TERMS'):
                    for goterm in obj.GO_TERMS.split('; '):
                        goterm = goterm.replace('|', '')
                        xref = XRef(goterm)
                        self.xrefs.append(xref)
                        self.edges.append(CreateEdge(peptide, xref, 'EVIDENCE'))
                        db_obj = self.db_checker("GO")
                        self.edges.append(CreateEdge(xref, db_obj, 'LINK_TO'))

                # creating edge to gene encoding the peptide
                obj.links_to_gene(peptide, self)

                if hasattr(obj, 'MODIFIED_FORM'):
                    uid = obj.MODIFIED_FORM
                    mod_obj = datfile.data[uid]
                    modification = Unspecified(
                        name=mod_obj.attr_check("COMMON_NAME", uid), uid=uid)
                    self.other_nodes.append(modification)
                    self.name_to_terms(modification)

                    # creating edges peptide -[FORMS]-> modification
                    self.edges.append(
                        CreateEdge(peptide, modification, 'FORMS'))

            for uid in complexes:
                obj = datfile.data[uid]

                # checking if a complex already in self.complexes
                complex_check = [i for i in self.complexes if i.uid == uid]
                if len(complex_check) == 0:
                    complex_obj = Complex(uid=uid,
                                          name=obj.attr_check("COMMON_NAME", uid),
                                          molecular_weight_kd=obj.attr_check("MOLECULAR_WEIGHT_KD"))
                    self.complexes.append(complex_obj)
                    self.name_to_terms(complex_obj)
                else:
                    complex_obj = complex_check[0]

                # creating Terms for complex name synonyms
                obj.links_to_synonyms(complex_obj, self)

                # creating XRef and DB node for GO-terms
                if hasattr(obj, 'GO_TERMS'):
                    for goterm in obj.GO_TERMS.split('; '):
                        goterm = goterm.replace('|', '')
                        xref = XRef(goterm)
                        self.xrefs.append(xref)
                        self.edges.append(CreateEdge(complex_obj, xref, 'EVIDENCE'))
                        db_obj = self.db_checker("GO")
                        self.edges.append(CreateEdge(xref, db_obj, 'LINK_TO'))

                if hasattr(obj, "COMPONENTS"):
                    db_objects = self.rnas + self.compounds + \
                                 self.oligopeptides + self.polypeptides + \
                                 self.complexes

                    comps = obj.COMPONENTS

                    for comp in comps.split("; "):
                        components = [c for c in db_objects if c.uid == comp]

                        if len(components) == 1:
                            component = components[0]

                        # if a component is a complex, that hasn't been created yet...
                        elif len(components) == 0 and comp[:4] == 'CPLX':
                            component = Complex(uid=comp,
                                                name=obj.attr_check("COMMON_NAME", comp),
                                                molecular_weight_kd=obj.attr_check("MOLECULAR_WEIGHT_KD"))
                            self.complexes.append(component)
                            self.name_to_terms(component)

                        # if a component is something else
                        else:
                            print "The complex with %s is assigned before " \
                                  "its components (%s)!" % (complex_obj.uid, comp)
                            continue
                        self.edges.append(
                            CreateEdge(component, complex_obj, 'PART_OF'))


            print "A list with %d oligopeptides has been created!\n" \
                    "A list with %d polypeptides has been created!\n" \
                    "A list with %d protein-complexes has been " \
                    "created!" % (len(self.oligopeptides),
                                  len(self.polypeptides), len(self.complexes))

    def protein_features_dat(self):
        """
        The method creates nodes for protein features (sequence conflicts and
        other) and links to proteins/peptides.
        """
        datfile = self._read_dat('protein-features.dat')
        if datfile is not None:
            for uid in datfile.names:
                obj = datfile.data[uid]
                protfeature = ProtFeature(uid=uid, type=obj.TYPES,
                                          residues=obj.feature_location(),
                                          source=obj.attr_check("DATA_SOURCE"),
                                          comment=obj.attr_check("COMMENT"),
                                          alternate_sequence=obj.attr_check("ALTERNATE_SEQUENCE"))
                self.protfeatures.append(protfeature)

                # creating an edge to the protein
                # Protein -[:HAS_FEATURE]-> ProtFeature
                obj.links_to_protein(protfeature, self)
            print "A list with %d protein features has been " \
                  "created!" % len(self.protfeatures)

    def sigma_factors(self):
        """
        The method creates nodes for Sigma-factors.
        """
        datfile = self._read_dat('proteins.dat')
        if datfile is not None:
            i = 0
            for uid in datfile.names:
                obj = datfile.data[uid]
                if "Sigma-Factors" in obj.TYPES.split('; '):
                    if hasattr(obj, "SYNONYMS"):
                        name = [s for s in obj.SYNONYMS.split('; ')
                                if s[:3] == 'Sig']
                        try:
                            name = name[0]
                        except:
                            name = obj.SYNONYMS.split('; ')[-1]
                    else:
                        name = obj.COMMON_NAME
                    sigma = SigmaFactor(name=name)
                    self.proteins.append(sigma)
                    self.name_to_terms(sigma)

                    # creating an edge to a poplypetide/complex
                    protein = [p for p in (self.polypeptides + self.complexes)
                               if p.uid == uid]
                    if len(protein) != 0:
                        self.edges.append(
                            CreateEdge(protein[0], sigma, 'IS_A'))

                    # creating edges to promoters
                    if hasattr(obj, "RECOGNIZED_PROMOTERS"):
                        for uid in obj.RECOGNIZED_PROMOTERS.split('; '):
                            pro = [p for p in self.promoters if p.uid == uid]
                            if len(pro) == 0:
                                continue
                            self.edges.append(
                                CreateEdge(sigma, pro, 'RECOGNIZES'))

                    i += 1
            print "%d Sigma-factors have been created!" % i

    def enzrxns_dat(self):
        """
        The method creates nodes for enzymes.
        """
        datfile = self._read_dat('enzrxns.dat')
        if datfile is not None:
            i = 0
            for uid in datfile.names:
                obj = datfile.data[uid]
                enzyme = Enzyme(uid=uid,
                                name=obj.attr_check("COMMON_NAME", uid))
                self.proteins.append(enzyme)
                self.name_to_terms(enzyme)

                # creating edge to a polypeptide or complex
                protein = [p for p in (self.polypeptides + self.complexes)
                           if p.uid == obj.ENZYME]
                if len(protein) == 0:
                    continue
                self.edges.append(CreateEdge(protein[0], enzyme, 'IS_A'))
                i += 1
            print "%d enzymes have been created!" % i

    def protcplxs_col(self):
        """
        The method assignes subunit composition information to complexes.
        """
        try:
            f = file(self.path + "protcplxs.col", 'r')
            data = f.readlines()
            f.close()
            for line in data:
                if line[0] == '#' or line[:3] == "UNI":
                    continue
                chunks = line.replace('\n', '').split('\t')
                complex_obj = [c for c in self.complexes
                               if c.uid == chunks[0]]
                if len(complex_obj) == 0 or len(chunks[-1]) == 0:
                    continue
                setattr(complex_obj[0], "subunit_composition", chunks[-1])
        except:
            print "There is no protcplxs.col file! Let's skip it..."

    def regulation_dat(self):
        """
        The method creates nodes for four subclasses of the RegulationEvent:
        Attenuation, EnzymeRegulation, TranscriptionRegulation,
        TranslationRegulation; and edges to Enzymes, Mechanisms, BSs,
        """
        datfile = self._read_dat('regulation.dat')
        if datfile is not None:
            att = ["Protein-Mediated-Attenuation",
                   "Small-Molecule-Mediated-Attenuation",
                   "RNA-Mediated-Attenuation"]
            transl = ["Protein-Mediated-Translation-Regulation",
                      "RNA-Mediated-Translation-Regulation",
                      "Compound-Mediated-Translation-Regulation"]
            for uid in datfile.names:
                obj = datfile.data[uid]
                if obj.TYPES == "Regulation-of-Enzyme-Activity":
                    reg = EnzymeRegulation(uid=uid,
                                           comment=obj.attr_check("COMMENT"),
                                           ki=obj.attr_check("KI"))
                elif obj.TYPES == "Transcription-Factor-Binding":
                    reg = TranscriptionRegulation(uid=uid,
                                                  comment=obj.attr_check("COMMENT"))
                elif obj.TYPES in att:
                    if hasattr(obj, "ANTITERMINATOR-END-POS") and \
                            hasattr(obj, "ANTITERMINATOR-START-POS"):
                        pos1 = [int(obj.ANTITERMINATOR_START_POS),
                                int(obj.ANTITERMINATOR_END_POS)]
                    else:
                        pos1 = None
                    if hasattr(obj, "ANTI-ANTITERM-START-POS") \
                            and hasattr(obj, "ANTI-ANTITERM-END-POS"):
                        pos2 = [int(obj.ANTI_ANTITERM_START_POS),
                                int(obj.ANTI_ANTITERM_END_POS)]
                    else:
                        pos2 = None
                    reg = Attenuation(uid=uid,
                                      comment=obj.attr_check("COMMENT"),
                                      antiterminator_pos=pos1,
                                      antiantiterminator_pos=pos2)
                elif obj.TYPES in transl:
                    reg = TranslationRegulation(uid=uid,
                                                comment=obj.attr_check("COMMENT"))
                else:
                    warnings.warn("Unexpected regulation types! "
                                  "Let's skip it... \n"
                                  "The object uid %s" % uid)
                    continue
                self.regulation_events.append(reg)

                # creating edges to objects participating in the regulation
                # (RegulationEvent)-[:REPRESSES/ACTIVATES/UNKNOWN]->(Node)
                obj.links_to_regulated_entity(reg, self)

                # creating edges to regulators
                # (Node)-[:PARTICIPATES_IN]->(RegulationEvent)
                obj.links_to_regulator(reg, self)

            print "A list with %d regulation events has been " \
                  "created!" % len(self.regulation_events)

    def create_compartments(self):
        """
        The method creates nodes for compartments.
        """
        common_names = ['Cytosol', 'Extracellular space', 'Cell wall',
                        'Periplasmic space', 'Cell envelope',
                        'Plasma membrane', 'Unknown']
        cco_names = ['CCO-CYTOSOL', 'CCO-EXTRACELLULAR', 'CCO-CW-BAC-POS',
                     'CCO-PERI-BAC', 'CCO-CE-BAC-POS', 'CCO-PM-BAC-POS',
                     'UNKNOWN']
        for cname, cco in zip(common_names, cco_names):
            self.compartments.append(Compartment(name=cname, uid=cco))

    def reactions_dat(self):
        """
        The method creates nodes for molecular reactions, reactants and edges
        to compartments.
        """
        datfile = self._read_dat('reactions.dat')
        if datfile is not None:
            other_nodes_num = len(self.other_nodes)
            for uid in datfile.names:
                obj = datfile.data[uid]

                # constructing the reaction formula
                if hasattr(obj, "LEFT") and hasattr(obj, "RIGHT"):
                    formula = obj.make_formula()
                    reaction = Reaction(uid=uid, type=obj.attr_check("TYPES"),
                                        formula=formula)
                    self.reactions.append(reaction)

                    # creating edges to enzymes catalyzing the reaction
                    # (Enzyme)-[:CATALYZES]->(Reaction)
                    obj.links_to_enzymes(reaction, self)

                    # creating edges to reaction name synonyms
                    obj.links_to_synonyms(reaction, self)

                    # searching for exact match of name/uid of objects
                    obj.links_to_reactants(reaction, self)

                    # creating edges to EC numbers
                    if hasattr(obj, 'EC_NUMBER'):
                        for ecnum in obj.EC_NUMBER.split('; '):
                            xref = XRef(ecnum)
                            self.xrefs.append(xref)
                            self.edges.append(CreateEdge(reaction, xref, 'EVIDENCE'))
                            db_obj = self.db_checker("ExPASy-ENZYME")
                            self.edges.append(CreateEdge(xref, db_obj, 'LINK_TO'))

                else:
                    continue
            print "A list with %d reactions has been created!\n" \
                  "There were %d unknown components, nodes were created " \
                  "for them..." % (len(self.reactions),
                                   len(self.other_nodes) - other_nodes_num)

    def pathways_dat(self):
        """
        The method creates nodes for pathways and edges to/from them.
        """
        datfile = self._read_dat('pathways.dat')
        if datfile is not None:
            superpath = []
            for uid in datfile.names:
                obj = datfile.data[uid]
                if obj.TYPES == 'Super-Pathways' or \
                                'Super-Pathways' in obj.TYPES.split('; '):
                    superpath.append(uid)
                    continue
                pathway = Pathway(uid=uid,
                                  name=obj.attr_check("COMMON_NAME", uid),
                                  type=obj.TYPES,
                                  reaction_layout=obj.attr_check("REACTION_LAYOUT"))
                self.pathways.append(pathway)
                self.name_to_terms(pathway)                
               
                # creating edges to reactions in the pathway
                # (Reaction) - [:PART_OF] -> (Pathway)
                obj.links_to_reactions(pathway, self)

                # creating links between reactions based on "PREDECESSORS"
                # attribute (Reaction) -[:NEXT]-> (Reaction)
                obj.reactions_order(self)

                # creating Terms for pathway name synonyms
                obj.links_to_synonyms(pathway, self)

            for uid in superpath:
                obj = datfile.data[uid]
                pathway = Pathway(uid=uid,
                                  name=obj.attr_check("COMMON_NAME", uid),
                                  type=obj.TYPES,
                                  reaction_layout=obj.attr_check("REACTION_LAYOUT"))
                self.pathways.append(pathway)
                self.name_to_terms(pathway) 

                # creating edges to reactions in the pathway
                # (Reaction) - [:PART_OF] -> (Pathway)
                # (Pathway) - [:PART_OF] -> (Pathway)
                obj.links_to_reactions(pathway, self)

                # creating Terms for pathway name synonyms
                obj.links_to_synonyms(pathway, self)
                
            print "A list with %d pathways has been " \
                  "created!" % len(self.pathways)

    def pathways_col(self):
        """
        The method creates nodes for genes acting in pathways.
        """
        if len(self.genes) == 0:
            print "There is no information about genes in the MetaCyc " \
                  "object! Let's skip it... \n"
        elif len(self.pathways) == 0:
            print "There is no information about pathways in the MetaCyc " \
                  "object! Let's skip it... \n"
        else:
            try:
                f = file(self.path + "pathways.col", 'r')
                data = f.readlines()
                f.close()
                if len(self.pathways) != 0:
                    for line in data:
                        if line[0] == '#' or line[:3] == "UNI":
                            continue
                        chunks = line.replace('\n', '').split('\t')
                        genes_uids = [chunk for chunk in chunks[2:]
                                      if chunk != '']
                        genes = [g for g in self.genes if g.uid in genes_uids]
                        pathway = [p for p in self.pathways
                                   if p.uid == chunks[0]][0]

                        # creating edges (Gene) -[ACTS_IN]-> (Pathway)
                        for gene in genes:
                            self.edges.append(
                                CreateEdge(gene, pathway, 'ACTS_IN'))

            except:
                print "There is no pathways.col file! Let's skip it..."

    def protseq_fsa(self):
        """
        The method adds sequences to peptide objects.
        """
        handle = open(self.path + "protseq.fsa", "rU")
        for record in SeqIO.parse(handle, "fasta"):
            uid = record.id.split("|")[2]
            protein = [p for p in (self.polypeptides + self.oligopeptides)
                       if p.uid == uid]
            if len(protein) == 0:
                continue
            else:
                protein[0].seq = str(record.seq)
        handle.close()

    def transporters_col(self):
        """
        The method creates nodes for proteins-transporters.
        """
        try:
            f = file(self.path + "transporters.col", 'r')
            data = f.readlines()
            f.close()
            for line in data:
                if line[0] == '#' or line[:3] == "UNI":
                    continue
                chunks = line.replace('\n', '').split('\t')
                transporter = Transporter(name=chunks[1], reaction=chunks[2])
                self.proteins.append(transporter)
                self.name_to_terms(transporter)
                protein = [p for p in (self.oligopeptides + self.polypeptides +
                                       self.complexes) if p.uid == chunks[0]]
                if len(protein) == 0:
                    continue

                # creating edges (Peptide) -[:IS_A]-> (Transporter)
                # creating edges (Complex) -[:IS_A]-> (Transporter)
                self.edges.append(CreateEdge(protein[0], transporter, 'IS_A'))
            transnum = len([p for p in self.proteins
                            if p.__class__.__name__ == "Transporter"])
            print "%d protein-transporters have been created!" % transnum
        except:
            print "There is no pathways.col file! Let's skip it..."

    def summary(self, filename="summary.txt"):
        """
        The method creates a file with MetaCyc object content overview.
        It also makes a number of tests to check graph characteristics
        and values formats.
        """
        # total number of nodes
        nodesnum = len(self.genes) + len(self.xrefs) + len(self.dbs)\
                   + len(self.terms) + len(self.rnas) + len(self.terminators)\
                   + len(self.promoters) + len(self.BSs) + len(self.TUs)\
                   + len(self.compounds) + len(self.polypeptides)\
                   + len(self.oligopeptides) + len(self.proteins)\
                   + len(self.complexes) + len(self.protfeatures)\
                   + len(self.regulation_events) + len(self.reactions)\
                   + len(self.pathways) + len(self.other_nodes)

        # protein groups numbers
        len_sf = len([p for p in self.proteins
                      if p.__class__.__name__ == "SigmaFactor"])
        len_e = len([p for p in self.proteins
                     if p.__class__.__name__ == "Enzyme"])
        len_t = len([p for p in self.proteins
                     if p.__class__.__name__ == "Transporter"])

        # file header
        result = "OrgID: %s \tOrganism: %s \tVersion: %s \tRelease date: %s\n" \
                 "\n" \
                 "Total number of nodes:\t %d\n" \
                 "Total number of edges:\t %d\n" \
                 "\n" % (self.orgid, self.organism, self.version, self.release,
                         nodesnum, len(self.edges))

        # constructing tables
        # Elements of the Genome
        table1 = tabulate([[len(self.genes), len(self.rnas),
                            len(self.promoters), len(self.terminators),
                            len(self.BSs), len(self.TUs)]],
                          headers=("Genes", "RNAs", "Promoters", "Terminators",
                                   "DNA Binding sites", "TUs"),
                          tablefmt="grid")

        # Compounds and various Proteins
        table2 = tabulate([[len(self.compounds), len(self.oligopeptides),
                            len(self.polypeptides), len_sf, len_e, len_t,
                            len(self.complexes)]],
                          headers=("Compounds", "Oligopeptides",
                                   "Polypeptides", "Sigma-Factors", "Enzymes",
                                   "Transporters", "Complexes"),
                          tablefmt="grid")

        # Regulation and Pathways
        table3 = tabulate([[len(self.regulation_events), len(self.reactions),
                            len(self.pathways), len(self.reactants),
                            len(self.compartments)]],
                          headers=("Regulation events", "Reactions",
                                   "Pathways", "Reactants", "Compartments"),
                          tablefmt="grid")

        # Special nodes
        table4 = tabulate([[len(self.terms), len(self.dbs), len(self.xrefs),
                            len(self.other_nodes)]],
                          headers=("Terms", "DBs", "XRefs", "Unspecified"),
                          tablefmt="grid")

        # WRITTING
        f = file(filename, "w")

        # information about release and number of elements: nodes and edges
        f.write(result)
        f.write("\nElements of the Genome\n")
        f.write(table1)
        f.write("\n\nCompounds and various Proteins\n")
        f.write(table2)
        f.write("\n\nRegulation and Pathways\n")
        f.write(table3)
        f.write("\n\nSpecial nodes and other\n")
        f.write(table4)

        # Results of a few tests for content
        f.write("\n\n\nContent tests\n"
                "---------------------\n\n")
        f.write(_DuplicateTest(self).result)
        f.write(_LonelyNodesTest(self).result)
        f.write(_NegativeCoordinatesTest(self).result)
        f.write(_StrandStringCheck(self).result)
        f.write(_XRefIDCheck(self).result)

        f.close()
        print "The txt-file with results has been created!"

    def make_graph(self, filename = 'metacyc'):
        """
        The method constructs a networkX graph structure from a MetaCyc object.
        All duplicates are deleting during graph constructuin.
        It might be very slow for objects with a great number of edges.
        SHOULD BE REWRITTEN!
        """
        j = 0
        allnodes = self.genes + self.xrefs + self.dbs + self.terms + \
                   self.rnas + self.terminators + self.promoters + \
                   self.BSs + self.TUs + self.compounds + \
                   self.polypeptides + self.oligopeptides + \
                   self.proteins + self.complexes + \
                   self.protfeatures + self.regulation_events + \
                   self.reactions + self.pathways + self.other_nodes + \
                   self.reactants + self.compartments

        graph = nx.DiGraph()

        # creating nodes
        nodes_dict = {}
        for node, i in zip(allnodes, xrange(0, len(allnodes))):
            mydict = node.__dict__
            for key in mydict.keys():
                if mydict[key] is None:
                    del mydict[key]
            graph_nodes = [g[1] for g in graph.nodes(data=True)]
            if mydict in graph_nodes:
                nodes_dict[str(i)] = graph_nodes.index(mydict)
            else:
                graph.add_node(j, mydict)
                nodes_dict[str(i)] = j
                j += 1
        print "Nodes done!"

        # creating edges
        i = 0
        edges = list(set(self.edges))
        for edge in edges:
            i += 1
            if i % 10000 == 0:
                print i
            try:
                i_source = allnodes.index(edge.source)
                i_target = allnodes.index(edge.target)
                graph.add_edge(nodes_dict[str(i_source)], nodes_dict[str(i_target)], label=edge.label)
            except:
                warnings.warn("Can't find nodes for an edge! Let's skip it!")
                continue

        nx.write_graphml(graph, filename + '.graphml')
        nx.write_gml(graph, filename + '.gml')
        return graph

    def organism_and_links(self, name, genome_status, parts):
        organism = Organism(name=name, genome_status=genome_status, parts=parts)
        self.other_nodes.append(organism)
        nodes = self.genes + self.rnas + self.promoters + self.terminators + \
                self.oligopeptides + self.polypeptides + self.proteins + \
                self.BSs + self.TUs + self.complexes + self.compartments
        for node in nodes:
            self.edges.append(CreateEdge(node, organism, 'PART_OF'))

###############################################################################


class _Test():
    """
    The class for various tests.
    """
    def __init__(self):
        self.result = None
        self.passed = None

###############################################################################


class _DuplicateTest(_Test):
    """
    The class is specially created to check content of a Metacyc object.
    Are there duplicated nodes or edges? The test will be passed only
    if there are no duplicates at all.
    """
    def __init__(self, metacyc):
        _Test.__init__(self)

        # checking nodes for duplicates
        allnodes = metacyc.genes + metacyc.xrefs + metacyc.dbs + \
                   metacyc.terms + metacyc.rnas + metacyc.terminators + \
                   metacyc.promoters + metacyc.BSs + metacyc.TUs + \
                   metacyc.compounds + metacyc.polypeptides + \
                   metacyc.oligopeptides + metacyc.proteins + \
                   metacyc.complexes + metacyc.protfeatures + \
                   metacyc.regulation_events + metacyc.reactions + \
                   metacyc.pathways + metacyc.other_nodes

        duplicates = len(allnodes) - len(set(allnodes))

        if duplicates == 0:
            self.passed = "Yes"
            self.result = 'Duplicate test for nodes:\t passed\n'
        else:
            self.passed = "No"
            self.result = 'Duplicate test for nodes:\t not passed\n' \
                          'There were about %d duplicates!\n' % duplicates

        # checking edges for duplicates
        duplicates = len(metacyc.edges) - len(set(metacyc.edges))

        if duplicates == 0 and self.passed == "Yes":
            self.passed = "Yes"
            self.result = self.result + 'Duplicate test for edges:\t ' \
                                        'passed\n\n'
        elif duplicates == 0 and self.passed == "No":
            self.result = self.result + 'Duplicate test for edges:\t ' \
                                        'passed\n\n'
        else:
            self.passed = "No"
            self.result += 'Duplicate test for edges:\t not passed\n' \
                          'There were about %d duplicates!\n\n' % duplicates

###############################################################################


class _LonelyNodesTest(_Test):
    """
    The class is specially created to check content of a Metacyc object.
    Are there lonely nodes without connection? The test will be passed
    only if there are no lonely nodes.
    """
    def __init__(self, metacyc):
        _Test.__init__(self)
        targets = [e.target for e in metacyc.edges]
        sources = [e.source for e in metacyc.edges]

        allnodes = metacyc.genes + metacyc.xrefs + metacyc.dbs + \
                   metacyc.terms + metacyc.rnas + metacyc.terminators + \
                   metacyc.promoters + metacyc.BSs + metacyc.TUs + \
                   metacyc.compounds + metacyc.polypeptides + \
                   metacyc.oligopeptides + metacyc.proteins + \
                   metacyc.complexes + metacyc.protfeatures + \
                   metacyc.regulation_events + metacyc.reactions + \
                   metacyc.pathways + metacyc.other_nodes

        diff = set(allnodes) - set(targets + sources)

        if len(diff) == 0:
            self.passed = 'Yes'
            self.result = 'Testing lonely nodes:\t passed\n\n'
        else:
            self.passed = 'No'
            self.result = 'Testing lonely nodes:\t not passed\n' \
                          'There were %d nodes without any connection!\n' \
                          '\n' % len(diff)

###############################################################################


class _NegativeCoordinatesTest(_Test):
    """
    The class is specially created to check content of a Metacyc object.
    Are there Feature objects with negative coordinates? The test will be
    passed only if all location coordinates are positive.
    """
    def __init__(self, metacyc):
        _Test.__init__(self)

        features = metacyc.genes + metacyc.terminators + \
                   metacyc.promoters + metacyc.BSs

        neg_coord = ''
        i = 0
        for feature in features:
            if feature.start < 0 or feature.end < 0:
                neg_coord += 'uid: %s\t%s\tstart = %d ' \
                             'end = %d\n' % (feature.uid, feature.labels,
                                             feature.start, feature.end)
                i += 1

        if len(neg_coord) == 0:
            self.passed = 'Yes'
            self.result = 'Testing Features for negative coordinates:\t ' \
                          'passed\n\n'
        else:
            self.passed = 'No'
            self.result = 'Testing Features for negative coordinates:\t ' \
                          'not passed\n' \
                          'There were %d nodes with negative coordinates.\n' \
                          '%s\n\n' % (i, neg_coord)

###############################################################################


class _StrandStringCheck(_Test):
    """
    The class is specially created to check content of a Metacyc object.
    Are there Feature objects with unexpected strand names? The test will be
    passed only if strand names are forward/reverse, unknown, both, +/-,
    sense/antisense.
    """
    def __init__(self, metacyc):
        _Test.__init__(self)

        features = metacyc.genes + metacyc.terminators + \
                   metacyc.promoters + metacyc.BSs
        possible_names = ['forward', 'reverse', 'unknown', 'both', '+', '-',
                          'sense', 'antisense']
        wrong_strand = ''
        i = 0
        for feature in features:
            if feature.strand not in possible_names:
                wrong_strand += 'uid: %s\t%s\t' \
                                'strand: %d\n' % (feature.uid, feature.labels,
                                                  feature.strand)
                i += 1

        if len(wrong_strand) == 0:
            self.passed = 'Yes'
            self.result = 'Testing Features for strand names:\t passed\n\n'
        else:
            self.passed = 'No'
            self.result = 'Testing Features for strand names:\t not passed\n' \
                          'There were %d nodes with unexpected strand ' \
                          'values.\n%s\n\n' % (i, wrong_strand)

###############################################################################


class _XRefIDCheck(_Test):
    """
    The class is specially created to check content of a Metacyc object.
    Are there XRefs objects with unexpected ids? The test will be
    passed only if all ids don't have spaces and symbols other then -.
    """
    def __init__(self, metacyc):
        _Test.__init__(self)

        wrong_id = ''
        i = 0
        for xref in metacyc.xrefs:
            if re.search(r'[^\w*\-.]', xref.id) is not None:
                wrong_id += 'id: %s\n' % xref.id
                i += 1

        if len(wrong_id) == 0:
            self.passed = 'Yes'
            self.result = 'Testing XRefs for id format:\t passed\n\n'
        else:
            self.passed = 'No'
            self.result = 'Testing XRefs for id format:\t not passed\n' \
                          'There were %d XRefs with unexpected id strings.\n' \
                          '%s\n\n' % (i, wrong_id)

###############################################################################


class _GenesTest(_Test):
    """
    The class is specially created to check information in a list of genes.
    Do genes have the same strand?
    Are they located closely to each other?
    """
    def __init__(self, geneslist):
        _Test.__init__(self)
        if not isinstance(geneslist, list):
            raise TypeError("The geneslist must be a list!")
        if len(geneslist) == 0:
            raise ValueError("The geneslist is empty!")
        types = list(set([g.__class__.__name__ for g in geneslist]))
        if types != ["Gene"]:
            raise TypeError("The geneslist must contain only genes!")
        self.geneslist = geneslist
        self.strand = "unknown"

        # If there is only one gene, we assign its strand information
        if len(geneslist) == 1:
            self.strand = geneslist[0].strand
            self.passed = "Yes"

        # For a number of genes -  we test them
        else:
            # Testing strand: is it the same or not?
            strand = list(set([g.strand for g in geneslist]))
            if len(strand) == 1:
                self.strand = strand[0]

            # Testing location: the difference between neighbouring genes
            # is less than 3000 bp
            starts = sorted([g.start for g in geneslist])[::-1]
            diff = [i - j for i, j, in zip(starts[:-1], starts[1:])]
            if max(diff) > 3000:
                self.passed = "No"
            else:
                self.passed = "Yes"

###############################################################################

#import doctest
#doctest.testfile("metacyc_tests_pathway_col.txt")
