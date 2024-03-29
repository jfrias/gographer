from xml.sax.handler import ContentHandler
from xml.sax import make_parser
from networkx import DiGraph
from GONode import *

class GOOboXmlHandler (ContentHandler):
    """ This class handles parsing the OBO XML file from the Gene Ontology Consortium and populate
        a GOGraph object passed to the constructor of this class.  The handler only use the nodes from
        the namespace (a branch of the GO) specificed in the GOGraph object.  Currently, only the edges
        that reflect IS_A relationship between a pair of GO terms are parsed and added into the GOGraph
    """
    ## Constructor.
    # @param    goGraph A reference to a GOGraph object to which the parsed GO nodes and edges
    #                   will be added.
    def __init__(self, goGraph):
        self.namespace = goGraph.getNameSpace()
        self.graph = goGraph
        self.inTypedef = False
        self.goNode = GONode()
        self.obsolete = False

    ## Create a new element
    # @param	name Set to either "term" if a new node is to be created or "typedef" if declaring a typedef is desired
	# @param	attributes
    def startElement(self, name, attributes):
        self.cdata = ""
        
        if name == "term":
            self.goid = None
            self.name = None
            self.description = None
            self.obsolete = False            
            self.elementNamespace = None
            self.parents = []
            self.goNode = GONode()
            
        elif name == "typedef":
            self.inTypedef = True
    ## Modify an element based on input
    # @param	name A string set to "id", "namespace", "is_a", "is_obsolete", "typedef", "term", "name", or "defstr"
	#				 depending on the desired changes
    def endElement(self, name):
        if name == "id":
            self.goid = self.cdata.strip()
            self.goNode.setGOID(self.goid)

        elif name == "namespace":
            self.elementNamespace = self.cdata.strip()
            self.goNode.setNamespace(self.elementNamespace)

        elif name == "is_a" and not self.inTypedef:
            self.parents.append(self.cdata.strip())
            self.goNode.setParents(self.parents)

        elif name == "is_obsolete" and self.cdata.strip() == "1":
            self.obsolete = True
            self.goNode.setObsolete(self.obsolete)

        elif name == "typedef":
            self.inTypedef = False

        elif name == "term" and not self.obsolete and self.namespace == self.elementNamespace:
            '''# not sure what the follow trying to do
            if not self.elementNamespace in self.namespaces.keys():
                self.namespaces[self.elementNamespace] = []

            self.namespaces[self.elementNamespace].append(self.goid)
            '''
            
            self.graph.add_node(self.goid, data=self.goNode)
            for parent in self.parents:
                '''
                currently the only edges being added are to the children of the current node
                '''
                if not parent in self.graph:
                    self.graph.add_node(parent, data=GONode(parent))
                self.graph.add_edge(parent, self.goid, relationship="parent_of")

        elif name == "name" and not self.obsolete and not self.goid == None and self.goNode.name == None: # and not self.name.has_key(self.goid):
            self.goNode.setName(self.cdata.strip())

        elif name == "defstr" and  not self.obsolete and not self.goid == None: # and not self.name.has_key(self.goid):
            self.goNode.setDescription(self.cdata.strip())

    ## Append additional data
	# @param	data Additional data to be appended
    def characters(self, data):
        self.cdata += data
