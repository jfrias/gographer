from xml.sax.handler import ContentHandler
from xml.sax import make_parser
from networkx import DiGraph
from GOOboXmlHandler import *
import cPickle

class GOGraph(DiGraph):
    ## Constructor.
    # @param    namespace     The branch of the GO ontology stored by this graph
    # @param    XMLFileName   The file contains the GO definition in
    #                         the format of the OBO in XML
    def __init__(self, namespace=None, GOOboXmlFileName=None):        
        DiGraph.__init__(self)

        self.namespace = namespace

        if GOOboXmlFileName != None:
            self.parseOboXml(GOOboXmlFileName)
        
    ## Parses the given OBO XML file and creates a GOGraph
    # @param    GOOboXMLFileName    The name of the OBO XML file to be parsed
    def parseOboXml(self, GOOboXmlFileName):
        parser = make_parser()
        handler = GOOboXmlHandler(self)
        parser.setContentHandler(handler)
        #try:
        f = open(GOOboXmlFileName, 'r')
        parser.parse(f)
        f.close()
        #except:
         #   print "Could not parse Obo XML file %s" % (GOOboXmlFileName)

    ## Returns the minimum depth of the node
    # @param    node    The node whose depth will be found and returned
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


    ## Returns the namespace of the graph
    def getNameSpace(self):
        return self.namespace;

    ## Save self using pickle protocol.
    # @param filename The filename of the pickle to use
    def savePickle(self, filename="gograph.pickle"):
        try:
            f = open(filename, 'wb')
            cPickle.dump(self, f, protocol = cPickle.HIGHEST_PROTOCOL)
            f.close()
        except:
            print "Could not pickle graph"

    ## Load a pickle from the filesystem
    # @param filename The location of the pickle file to load
    @classmethod
    def loadPickle (klass, filename="gograph.pickle"):
        try:
            f = open(filename, 'rb')
            g = cPickle.load(f)
            return g
        except:
            print "Could not load pickle"
            return
