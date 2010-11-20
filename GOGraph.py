from xml.sax.handler import ContentHandler
from xml.sax import make_parser
from networkx import DiGraph
from parseGOOboXml import *

class GOGraph(DiGraph):
    def __init__(self, aspect, XML):
        self.aspect = aspect
        self.parseOboXml(XML)

    def parseOboXml(self, XML):
        parser = make_parser()
        handler = parseGOOboXml()
        parser.setContentHandler(handler)
        try:
            f = open(XML, 'r')
            parser.parse(f)
            f.close()
            self.graph = handler.graph
        except:
            print "Could not parse OWL file %s" % (XML)

    def getLevel(self, node):
        parents = self.graph.predecessors(node)
        if len(parents) == 0:
            return 0
        
        found = False
        level = 1
        while not found:
            nextDepth = list()
            for i in parents:
                pred = self.graph.predecessors(i)
                if len(pred) == 0:
                    found = True
                    return level
                else:
                    nextDepth = list(set(nextDepth) | set(pred))
            parents = nextDepth
            level = level + 1


location = "./go_daily-termdb.obo-xml"

GOGraph = GOGraph("temp", location)

