## \package GOGeneGraph

from GOGraph import GOGraph
from networkx import topological_sort

class GOGeneGraph(GOGraph):
    ## Create a gene graph from a GOGraph
    # @param    gograph A GOGraph to base this graph off of
    # @param    assoc   The file containing gene association information
    def __init__(self, gograph, assoc=None):
        if not gograph.__class__.__name__ is 'GOGraph':
            msg =  "You did not give an instance of GOGraph to GOProteinGraph (it's a %s):" % gograph.__class__.__name__ \
                  + " any weights or directions will be ignored."
            
        GOGraph.__init__(self)
        self.add_edges_from(gograph.edges_iter(data=True))
        self.add_nodes_from(gograph.nodes_iter(data=True))
        self.geneToNode = dict()

        if assoc != None:
            self.parseAssocFile(assoc)
            self.propagateGenes()
            
    ## Parses the given association file and adds the gene information to the appropriate nodes
    # @param    assoc   The name of the association file to be parsed
    def parseAssocFile(self, assoc):
        types = ["protein"]
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

                    #Stores both the id and the qualifier as a tuple in the node's gene list
                    self.node[fields[4]]['data'].addGenes([(fields[1],fields[3])])

                    #Adds the gene to node association to the dictionary only if it isn't already in dictionary
                    if fields[1] not in self.geneToNode:
                        self.geneToNode[fields[1]] = set([fields[4]])
                    else:
                        self.geneToNode[fields[1]] = self.geneToNode[fields[1]].union([fields[4]])
            f.close
        except:
            print "Could not parse association file %s" % (assoc)
    
    ## Applies a weight to the graph
    def weight(self):
        raise NotImplementedError

    ## Returns a GOGenePubmedGraph version of itself
    # @param    gopubmedgraph   The GOPubmedGraph that this graph is to be merged with
    def toGOGenePubmedGraph(self):
        raise NotImplementedError

    ## Returns the associated genes for a given node
    # @param    goid    The GOID of the node to retrieve the associated genes from
    def getGenesByNode(self, goid):
        if goid in self:
            return self.node[goid]['data'].getGenes()
        else:
            print "Given GOID is not in the graph"
            raise

    ## Returns the associated propagated genes for a given node
    # @param    goid    The GOID of the node to retrieve the associated propagated genes from
    def getPropagatedGenesByNode(self, goid):
        if goid in self:
            return self.node[goid]['data'].getPropagatedGenes()
        else:
            print "Given GOID is not in the graph"
            raise


    ## Returns the associated nodes for a given gene
    # @param    geneid  The ID of the gene for which to retrieve the associated nodes
    def getNodesByGene(self, geneid):
        if geneid in self.geneToNode:
            return self.geneToNode[geneid]
        else:
            print "Given gene id is not in the graph"
            raise

    ## Propagate the genes associated with each node to its parents
    def propagateGenes(self):
        for node in self.nodes_iter():
            self.node[node]['data'].addPropagatedGenes(self.node[node]['data'].getGenes())

        sortedNodes = topological_sort(self)
        sortedNodes.reverse()

        # go through every nodes from bottom (end of sorted list)
        # propagate node.propPmids to parent nodes
        for node in sortedNodes:
            for ancestors in self.predecessors(node):
                self.node[ancestors]['data'].addPropagatedGenes(self.node[node]['data'].getPropagatedGenes())

    def removeGeneless(self):
        geneless = []
        for node in self:
            if len(self.getPropagatedGenesByNode(node)) == 0:
                geneless.append(node)
        self.remove_nodes_from(geneless)
