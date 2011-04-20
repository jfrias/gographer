from xml.sax.handler import ContentHandler
from xml.sax import make_parser
from networkx import DiGraph
from networkx import topological_sort
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
        try:
            f = open(GOOboXmlFileName, 'r')
            parser.parse(f)
            f.close()
        except:
            print "Could not parse Obo XML file %s" % (GOOboXmlFileName)

    ## Returns the minimum depth of the node
    # @param    goid    The GO ID of the node whose depth will be found and returned
    def getLevel(self, goid):
        parents = self.predecessors(goid)
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
    # @param    filename    The location of the pickle file to load
    @classmethod
    def loadPickle (klass, filename="gograph.pickle"):
        try:
            f = open(filename, 'rb')
            g = cPickle.load(f)
            return g
        except:
            print "Could not load pickle"
            return

    ## Get the description of a node
    # @param    goid   The GO ID of the node to get the description of
    def getNodeDescription(self, goid):
        if goid in self.nodes():
            return self.node[goid]['data'].getDescription()
        else:
            print 'Invalid goid'
            raise

    ## Set the description of a node
    # @param    goid   The GO ID of the node to get the description of
    # @param    descrip The description that will be assigned to the node
    def __setNodeDescription(self, goid, descrip):
        if goid in self.nodes():
            self.node[goid]['data'].setDescription(descrip)
        else:
            print 'Invalid goid'
            raise

    ## Get the description of a node
    # @param    goid   The GO ID of the node to get the description of
    def getNodeNamespace(self, goid):
        if goid in self.nodes():
            return self.node[goid]['data'].getNamespace()
        else:
            print 'Invalid goid'
            raise

    ## Set the description of a node
    # @param    goid   The GO ID of the node to get the description of
    # @param    namespace The namespace that will be assigned to the node
    def __setNodeNamespace(self, goid, namespace):
        if goid in self.nodes():
            self.node[goid]['namespace'] = namespace
        else:
            print 'Invalid goid'
            raise

    ## Calculates and stores the number of descendants each node in the graph has
    def calcDescendantCount(self):
        nodes = dict()
        sortedNodes = topological_sort(self)
        sortedNodes.reverse()
        for node in sortedNodes:
            if node not in nodes:
                ids = set()
            else:
                ids = nodes[node]

            self.node[node]['data'].descendantCount = len(ids)
            ids = ids.union([node])

            for parent in self.predecessors(node):
                if parent not in nodes:
                    nodes[parent] = ids
                else:
                    nodes[parent] = nodes[parent].union(ids)

    ## Return the number of descendants of the given node. Calculates the count if it had not already been done so
    # @param    goid    The GO ID of the node to return the descendant count of
    def getDescendantCount(self, goid):
        if goid in self:
            count = self.node[goid]['data'].getDescendantCount()
            if not count:
                self.calcDescendantCount()
                return self.node[goid]['data'].getDescendantCount()
            else:
                return count

