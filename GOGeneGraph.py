## \package GOGeneGraph

from GOGraph import GOGraph

class GOGeneGraph(GOGraph):
    ## Create a gene graph from a GOGraph
    # @param    gograph A GOGraph to base this graph off of
    # @param    assoc   The file containing gene association information
    def __init__(self, gograph, assoc):


    ## Applies a weight to the graph
    #
    def weight(self):

    ## Returns a GOGenePubmedGraph version of itself
    # @param    gopubmedgraph   The GOPubmedGraph that this graph is to be merged with
    def toGOGenePubmedGraph(self):
        
