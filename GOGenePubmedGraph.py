## @package GOGenePubmedGraph

from GOPubmedGraph import GOPubmedGraph
from GOGeneGraph import GOGeneGraph

class GOGenePubmedGraph(GOPubmedGraph, GOGeneGraph):
    ## Create a gene pubmed graph from a GOPubmedGraph and GOGeneGraph
    # @param    gopubmedgraph   The GOPubmedGraph that should be merged to form the gene pubmed graph
    # @param    gogenegraph The GOGeneGraph that should be merged to form the gene pubmed graph
    def __init__(self, gopubmedgraph, gogenegraph):
        raise NotImplementedError

    ## Weights the gene pubmed graph based off of the weight factor
    # @param    weightFactor    The factor by which the proteins and semantic distances should be weighted
    def weight(self, weightFactor):
        raise NotImplementedError
        
