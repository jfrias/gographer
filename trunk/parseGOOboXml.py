from xml.sax.handler import ContentHandler
from xml.sax import make_parser
from networkx import DiGraph
from GONode import *

class parseGOOboXml(ContentHandler):
    def __init__(self):
        self.aspects = {}
        self.relationships = {}
        self.inTypedef = False
        self.names = {}
        self.obsolete = False
        self.graph = DiGraph()
        self.goNode = GONode()
        
    def startElement(self, name, attributes):
        self.cdata = ""
        if name == "term":
            self.goid = None
            self.aspect = None
            self.parents = []
            self.obsolete = False
            self.goNode = GONode()
        elif name == "typedef":
            self.inTypedef = True
                
    def endElement(self, name):
        if name == "id":
            self.goid = self.cdata.strip()
            self.goNode.setGOID(self.goid)

        elif name == "namespace":
            self.aspect = self.cdata.strip()
            self.goNode.setNamespace(self.aspect)

        elif name == "is_a" and not self.inTypedef:
            self.parents.append(self.cdata.strip())
            self.goNode.setParents(self.parents)

        elif name == "is_obsolete" and self.cdata.strip() == "1":
            self.obsolete = True
            self.goNode.setObsolete(self.obsolete)

        elif name == "typedef":
            self.inTypedef = False

        elif name == "term" and not self.obsolete:
            self.relationships[self.goid] = self.parents

            if not self.aspect in self.aspects.keys():
                self.aspects[self.aspect] = []

            self.aspects[self.aspect].append(self.goid)

            self.graph.add_node(self.goid, data=self.goNode)
            for parent in self.relationships[self.goid]:
                self.graph.add_edge(parent, self.goid)

        elif name == "name" and not self.obsolete and not self.goid == None and not self.names.has_key(self.goid):
            self.names[self.goid] = self.cdata.strip()
            self.goNode.setName(self.names[self.goid])

    def characters(self, data):
        self.cdata += data
