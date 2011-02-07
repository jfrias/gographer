## \package GOGeneGraph

from GOGraph import GOGraph

class GOGeneGraph(GOGraph):
    ## Create a gene graph from a GOGraph
    # @param    gograph A GOGraph to base this graph off of
    # @param    assoc   The file containing gene association information
    def __init__(self, gograph, assoc):
        raise NotImplementedError

    ## Applies a weight to the graph
    #
    def weight(self):
        raise NotImplementedError

    ## Returns a GOGenePubmedGraph version of itself
    # @param    gopubmedgraph   The GOPubmedGraph that this graph is to be merged with
    def toGOGenePubmedGraph(self):
        raise NotImplementedError

    ## Returns the associated genes for a given node
    # @param    goid    The GOID of the node to retrieve the associated genes from
    def getGenesByNode(self, goid):
        raise NotImplementedError

    ## Returns the associated nodes for a given gene
    # @param    geneid  The ID of the gene for which to retrieve the associated nodes
    def getNodesByGene(self, geneid):
        raise NotImplementedError
