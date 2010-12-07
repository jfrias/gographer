from xml.sax.handler import ContentHandler
from xml.sax import make_parser
from networkx import DiGraph
from GOOboXmlHandler import *

class GOGraph(DiGraph):
    def __init__(self, namespace, GOOboXmlFileName):
        """ Constructor.
        Args:
            namespace      The branch of the GO ontology stored by this graph
            XMLFileName    The file contains the GO definition in the format
                           of the OBO in XML
        """
        DiGraph.__init__(self)
        
        self.namespace = namespace
        self.parseOboXml(GOOboXmlFileName)
        

    def parseOboXml(self, GOOboXmlFileName):
        parser = make_parser()
        handler = GOOboXmlHandler(self)
        parser.setContentHandler(handler)
        try:
            f = open(GOOboXmlFileName, 'r')
            parser.parse(f)
            f.close()
        except:
            print "Could not parse Obo XML file %s" % (XML)

    def getLevel(self, node):
        parents = self.predecessors(node)
        if len(parents) == 0:
            return 0
        
        found = False
        level = 1
        while not found:
            nextDepth = list()
            for i in parents:
                pred = self.predecessors(i)
                if len(pred) == 0:
                    found = True
                    return level
                else:
                    nextDepth = list(set(nextDepth) | set(pred))
            parents = nextDepth
            level = level + 1


    def getNameSpace(self):
        return self.namespace;
    
