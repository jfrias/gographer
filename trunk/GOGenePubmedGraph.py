## @package GOGenePubmedGraph

from GOPubmedGraph import GOPubmedGraph
from GOGeneGraph import GOGeneGraph
from networkx import topological_sort

class GOGenePubmedGraph(GOPubmedGraph, GOGeneGraph):
    ## Create a gene pubmed graph from a GOPubmedGraph and GOGeneGraph
    # @param    gopubmedgraph   The GOPubmedGraph that should be merged to form the gene pubmed graph
    # @param    gogenegraph The GOGeneGraph that should be merged to form the gene pubmed graph
    def __init__(self, gopubmedgraph, gogenegraph):
        GOGeneGraph.__init__(self, gogenegraph)
        #GOPubmedGraph.__init__(self, gopubmedgraph)
        self.add_edges_from(gopubmedgraph.edges_iter(data=True))
        self.add_edges_from(gogenegraph.edges_iter(data=True))
        self.add_nodes_from(gogenegraph.nodes_iter(data=True))
        for node in gopubmedgraph:
            wordVector = gopubmedgraph.node[node]['data'].getWordVector()
            self.node[node]['data'].wordVector = wordVector.copy()
            self.node[node]['data'].pmids = gopubmedgraph.node[node]['data'].pmids.copy()
            self.node[node]['data'].propPmids = gopubmedgraph.node[node]['data'].propPmids.copy()
        self.pubmedToNode = gopubmedgraph.pubmedToNode.copy()
        self.geneToNode = gogenegraph.geneToNode.copy()
        self.excludeEvidence = gopubmedgraph.excludeEvidence

    ## Weights the gene pubmed graph based off of the weight factor
    # @param    weighter    The weighing method by which the proteins and semantic distances should be weighted
    def weightNodes(self, weighter):
        sortedNodes = topological_sort(self)
        sortedNodes.reverse()

        for node in sortedNodes:
            for parent in self.predecessors(node):
                weight = weighter.makeWeighted(node, parent, self)
                self.edge[parent][node]['weight'] = weight
        