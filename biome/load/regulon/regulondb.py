from ...api import *

class RegulonDB():
    """

    """
    def __init__(self, ecoli_node, connection='http://localhost:7474/db/data/'):
        if not isinstance(ecoli_node, Node):
            raise TypeError('The ecoli_node argument must be a Node!')
        if not isinstance(connection, basestring):
            raise TypeError('The connection argument must be a string!')
        self.ecoli_node = ecoli_node
        self.connection = connection