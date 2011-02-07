## @package GOPubmedGraph

from GOGraph import GOGraph

class GOPubmedGraph(GOGraph):
    ## Create a pubmed graph from a GOGraph
    # @param    gograph A GOGraph to base this graph off of
    # @param    XMLFileName The file containing GO definitions and Pubmed ID
    #                       information in the format of the OBO in XML
    def __init__(self, gograph, XMLFileName):
        raise NotImplementedError

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
        raise NotImplementedError

    ## Returns the associated nodes for a given PubMed ID
    # @param    pubmed  The PubMed ID for which to retrieve the associated nodes:
    def getNodesByPubMed(self, pubmed):
        raise NotImplementedError

