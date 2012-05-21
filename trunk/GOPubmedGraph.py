## @package GOPubmedGraph

from GOGraph import GOGraph
from networkx import topological_sort
from utils import *
import math

class GOPubmedGraph(GOGraph):
    ## Create a pubmed graph from a GOGraph
    # @param    gograph A GOGraph to base this graph off of
    # @param    assoc   The file containing gene and pubmed association information
    # @param    excludeEvidence A list of the evidence codes that should be ignored
    def __init__(self, gograph, assoc=None, excludeEvidence =[], types = ["protein"]):
        if not gograph.__class__.__name__ is 'GOGraph':
            msg =  "You did not give an instance of GOGraph to GOProteinGraph (it's a %s):" % gograph.__class__.__name__ \
                  + " any weights or directions will be ignored."
            
        GOGraph.__init__(self)
        self.add_edges_from(gograph.edges_iter(data=True))
        self.add_nodes_from(gograph.nodes_iter(data=True))
        self.pubmedToNode = dict()
        self.excludeEvidence = excludeEvidence

        if assoc != None:
            self.parseAssocFile(assoc, types)

    ## Parses the given association file and adds the pubmed information to the appropriate nodes
    # @param    assoc   The name of the association file to be parsed
    def parseAssocFile(self, assoc, types = ["protein"]):
        nodes = self.nodes()
        try:
            f = open(assoc, 'r')
            for line in f:
                if line[0] == '!':
                    continue
                else:
                    fields = line.split('\t')
                    #Checks to make sure the go term is in the graph
                    if not fields[4] in nodes:
                        continue
                    #Checks to make sure the association is a type of interest
                    if not fields[11] in types:
                        continue

                    #Checks to make sure the evidence type isn't one that was listed to be excluded
                    if fields[6] in self.excludeEvidence:
                        continue

                    # Pmid's that start w/ P may bes provincial
                    references = filter(lambda x: x.startswith("PMID:") and not x.startswith("PMID:P"), fields[5].split('|'))

                    if(len(references) > 0):
                        for ref in references:
                            #Stores both the id and the qualifier in the node's pubmed list
                            self.node[fields[4]]['data'].addPMIDs([(ref[5:],fields[3])])

                            #Adds the pubmed to node association to the dictionary only if it isn't already in dictionary
                            if ref[5:] not in self.pubmedToNode:
                                self.pubmedToNode[ref[5:]] = set([fields[4]])
                            else:
                                self.pubmedToNode[ref[5:]] = self.pubmedToNode[ref[5:]].union([fields[4]])
            f.close
            self.propagatePMIDs()
        except:
            print "Could not parse association file %s" % (assoc)

    ## Applies a weight to the graph
    def weight(self):
        raise NotImplementedError

    ## Returns a GOGenePubmedGraph version of itself
    # @param    gogenegraph The GOGeneGraph that this graph is to be merged with
    def toGOGenePubmedGraph(self):
        raise NotImplementedError

    ## Returns the associated PubMed IDs for a given node
    # @param    goid    The GOID of the node to retrieve the associated PubMed IDs from
    def getPubMedByNode(self, goid):
        if goid in self:
            return self.node[goid]['data'].getPMIDs()
        else:
            print "Given GOID is not in the graph"
            raise

    ## Returns the associated propagated PubMed IDs for a given node
    # @param    goid    The GOID of the node to retrieve the associated propagated PubMed IDs from
    def getPropagatedPubMedByNode(self, goid):
        if goid in self:
            return self.node[goid]['data'].getPropagatedPMIDs()
        else:
            print "Given GOID is not in the graph"
            raise

    ## Returns the associated nodes for a given PubMed ID
    # @param    pubmed  The PubMed ID for which to retrieve the associated nodes:
    def getNodesByPubMed(self, pubmed):
        if pubmed in self.pubmedToNode:
            return self.pubmedToNode[pubmed]
        else:
            print "Given pubmed id is not in the graph"
            raise

    ## Propagate the PMIDs associated with each node to its parents
    def propagatePMIDs(self):
        for node in self.nodes_iter():
            self.node[node]['data'].addPropagatedPMIDs(self.node[node]['data'].getPMIDs())
        
        sortedNodes = topological_sort(self)
        sortedNodes.reverse()

        # go through every nodes from bottom (end of sorted list)
        # propagate node.propPmids to parent nodes
        for node in sortedNodes:
            for ancestors in self.predecessors(node):
                self.node[ancestors]['data'].addPropagatedPMIDs(self.node[node]['data'].getPropagatedPMIDs())

    ##Calculates the word vectors for all of the nodes in the graph
    # @param    corpus  The corpus that contains the information on the PubMed article
    # @param    tokenizer   The tokenizer function that will be used on the text, a simple tokenizer is used if none is given.
    #                       Should take a string as an input, and outputs a string with words that are lower case and separated by a space
    # @param    stemmer The stemmer function that will be used to stem the words, the porter stemmer is used if none is given
    #                   Takes a word as an input and reports a stemmed word as an output.
    # @param    stopwords   A list of stop words that will not be included in the word vector, either as a list or a StopwordList.
    #                       An empty list is used if no stop word list is given.
    def calculateWordVectors(self, corpus, tokenizer=Tokenizer().tokenize_word,
                            stemmer=PorterStemmer().stem, stopwords=[]):
        for node in self.nodes_iter():
            self.node[node]['data'].calculateWordVector(corpus, tokenizer, stemmer, stopwords)

    ##Remove nodes that don't have at least one associated PMID
    def removePMIDless(self):
        pmidless = []
        for node in self:
            if len(self.getPropagatedPubMedByNode(node)) == 0:
                pmidless.append(node)
        self.remove_nodes_from(pmidless)

    ##Temporary function for finding common ancestors
    '''
    def findCommonAncestors(self, terms):
        temp = DiGraph()
        temp.add_edges_from(self.edges(data=True))
        reverse = temp.reverse()
        
        ancestors = dict()
        for term in terms:
            for ancestor in shortest_path_length(reverse, term, weighted=True):
                if ancestor in ancestors:
                    ancestors[ancestor] = ancestors[ancestor].union([term])
                else:
                    ancestors[ancestor] = set([term])

        key = ancestors.keys()
        for ancestor in key:
            if len(ancestors[ancestor]) == 1:
		ancestors.pop(ancestor)

        return ancestors
    '''
