## @package GOGenePubmedGraph

from GOPubmedGraph import GOPubmedGraph
from GOGeneGraph import GOGeneGraph
from networkx import topological_sort
from networkx import Graph
from networkx import DiGraph
from heapq import *
from utils import keepGenes
from calcModel import *
from mergeGraph import *

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

    ## Parses the given association file and adds the gene and pubmed information to the appropriate nodes
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

                    if fields[6] in self.excludeEvidence:
                        continue

                    #Stores both the id and the qualifier as a tuple in the node's gene list
                    self.node[fields[4]]['data'].addGenes([(fields[1],fields[3])])

                    #Adds the gene to node association to the dictionary only if it isn't already in dictionary
                    if fields[1] not in self.geneToNode:
                        self.geneToNode[fields[1]] = set([fields[4]])
                    else:
                        self.geneToNode[fields[1]] = self.geneToNode[fields[1]].union([fields[4]])
                    
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
            self.propagateGenes()
        except:
            print "Could not parse association file %s" % (assoc)

    ## Weights the gene pubmed graph based off of the weight factor
    # @param    weighter    The weighing method by which the proteins and semantic distances should be weighted
    def weightNodes(self, weighter):
        sortedNodes = topological_sort(self)
        sortedNodes.reverse()

        for node in sortedNodes:
            for parent in self.predecessors(node):
                weight = weighter.makeWeighted(node, parent, self)
                self.edge[parent][node]['weight'] = weight

    ## Merges the graph using the selected method
    # @param    genes   List of genes that should be kept in the graph
    # @param    mode    Merging method to be used, default is 'IB'
    # @param    model   Model used for calculating probability
    # @param    maxProb The maximum allowed probability for a merging to occur, the default is 0.05
    # @param    maxMergedGeneCount  The maximum allowed number of merged genes for any given node, merging stops at that node if the max is reached. The default value is 200
    # @param    maxGeneAutoMerge    The maximum number of genes associated with a node for it to be automatically merged
    def makeMerged(self, genes, model, mode='IB', maxProb=0.05, maxMergedGeneCount=200, minGeneAutoMerge=5):
        keepGenes(self, genes)
        self.removeGeneless()
        self.removePMIDless()
        if mode == "IB":
            self.mergeGraphIB(model, maxProb, maxMergedGeneCount, minGeneAutoMerge)
        if mode == 'Mult':
            self.mergeGraphMult(model, maxProb, maxMergedGeneCount, minGeneAutoMerge)
        if mode == 'Prob':
            self.mergeGraphProb(model, maxProb, maxMergedGeneCount, minGeneAutoMerge)

    ## Merges the graph, and checks before each merge to ensure the probability of associated information loss is less than the set max probability
    # @param graph A weighted GOGenePubmedGraph that will be merged
    # @param model The model from which to calculate the probability
    # @param maxProb The max allowed probaiblity for a merging to occur, the default value is 0.05
    # @param maxMergedGeneCount The max allowed number of merged genes for any given node, merging stops at that node if the max is reached. The default value is 200
    def mergeGraphIB(self, model, maxProb=0.05, maxMergedGeneCount=200, minGeneAutoMerge=5):
        queue = []
        leafs = set()
        for node in self.nodes():
            if len(self.edges(node)) == 0:
                leafs.add(node)
                for parent in self.predecessors(node):
                    heappush(queue, (self.edge[parent][node]['weight'], (parent, node)))

        while len(queue) > 0:
            edge = heappop(queue)
            #Calculates information loss
            loss = self.node[edge[1][0]]['data'].getInfoLoss() + edge[0] + self.node[edge[1][1]]['data'].getInfoLoss()
            mergedNodes = self.node[edge[1][0]]['data'].getMergedCount() + self.node[edge[1][1]]['data'].getMergedCount() + 1
            mergedGenes = len(self.node[edge[1][0]]['data'].getMergedGenes().union(self.node[edge[1][1]]['data'].getMergedGenes()).union(self.getGenesByNode(edge[1][1])).union(self.getGenesByNode(edge[1][0])))
            genes = len(self.node[edge[1][1]]['data'].getMergedGenes().union(self.getGenesByNode(edge[1][1])))

            #Obtains probability
            if mergedGenes <= maxMergedGeneCount:
                prob = calcProb(loss, mergedGenes, model)
            else:
                prob = 1 + maxProb

            #Merges if probability is less than maxProb
            if prob < maxProb or genes <= minGeneAutoMerge:
                predecessors = self.predecessors(edge[1][1])
                self = mergeEdge(self, (edge[1][0], edge[1][1]))
                self.remove_node(edge[1][1])
                leafs.remove(edge[1][1])

                for pred in predecessors:
                    if len(self.edges(pred)) == 0:
                        self.removeEmptyNode(pred)

                queue = []
                for node in self.nodes():
                    if len(self.edges(node)) == 0:
                        leafs.add(node)
                        for parent in self.predecessors(node):
                            heappush(queue, (self.edge[parent][node]['weight'], (parent, node)))
        return self, leafs

	## Multi-Graph Merge
	# @param	model
	# @param	maxProb Always set to 0.05
	# @param	maxMergedGeneCount Always set to 200
	# @param	minGeneAutoMerge Always set to 5
    def mergeGraphMult(self, model, maxProb=0.05, maxMergedGeneCount=200, minGeneAutoMerge=5):
        queue = []
        leafs = set()
        for node in self.nodes():
            if len(self.edges(node)) == 0:
                leafs.add(node)
                for parent in self.predecessors(node):
                    geneCount = len(self.getGenesByNode(node))
                    if geneCount == 0:
                        geneCount = 1
                    heappush(queue, (self.edge[parent][node]['weight']*(geneCount+len(self.node[node]['data'].getMergedGenes())), (parent, node)))

        while len(queue) > 0:
            edge = heappop(queue)
            #Calculates information loss
            loss = self.node[edge[1][0]]['data'].getInfoLoss() + self.edge[edge[1][0]][edge[1][1]]['weight'] + self.node[edge[1][1]]['data'].getInfoLoss()
            mergedNodes = self.node[edge[1][0]]['data'].getMergedCount() + 1 + self.node[edge[1][1]]['data'].getMergedCount() + 1
            mergedGenes = len(self.node[edge[1][0]]['data'].getMergedGenes().union(self.node[edge[1][1]]['data'].getMergedGenes()).union(self.getGenesByNode(edge[1][1])).union(self.getGenesByNode(edge[1][0])))
            genes = len(self.node[edge[1][1]]['data'].getMergedGenes().union(self.getGenesByNode(edge[1][1])))

            #Obtains probability
            if mergedGenes <= maxMergedGeneCount:
                prob = calcProb(loss, mergedGenes, model)
            else:
                prob = 1 + maxProb

            #Merges if probability is less than maxProb
            if prob < maxProb or genes <= minGeneAutoMerge:
                predecessors = self.predecessors(edge[1][1])
                self = mergeEdge(self, (edge[1][0], edge[1][1]))
                self.remove_node(edge[1][1])
                leafs.remove(edge[1][1])

                for pred in predecessors:
                    if len(self.edges(pred)) == 0:
                        self.removeEmptyNode(pred)

                queue = []
                for node in self.nodes():
                    if len(self.edges(node)) == 0:
                        leafs.add(node)
                        for parent in self.predecessors(node):
                            geneCount = len(self.getGenesByNode(node))
                            if geneCount == 0:
                                geneCount = 1
                            heappush(queue, (self.edge[parent][node]['weight']*(geneCount+len(self.node[node]['data'].getMergedGenes())), (parent, node)))
                
        return self, leafs
	
	## Augmented Graph Merge
	# @param	model
	# @param	maxProb Always set to 0.05
	# @param	maxMergedGeneCount Always set to 200
	# @param	minGeneAutoMerge Always set to 5
	# @param	minLevel Always set to 0
    def mergeAugmented(self, model, maxProb=0.05, maxMergedGeneCount=200, minGeneAutoMerge=5, minLevel=0):
        #get undirected graph
        undirected = self.to_undirected()

        #calculate descendant dictionary
        descendantDict = dict()
        sortedNodes = topological_sort(self)
        sortedNodes.reverse()
        for node in sortedNodes:
            descendants = descendantDict.get(node, set())
            descendants.add(node)
            descendantDict[node] = descendants

            for parent in self.predecessors(node):
                if parent not in descendantDict:
                    descendantDict[parent] = set(descendants)
                else:
                    descendantDict[parent].update(descendants)

        #calculate 5th percentile weight
        sort = list()
        for edge in self.edges(data=True):
            sort.append(edge[2]['weight'])
        sort.sort()
        fifth = sort[int(math.ceil(len(sort))*0.05-1)]

        genes = self.geneToNode.keys()

        geneTuples = set()
        for gene in genes:
            geneTuples.add((gene,''))
            geneTuples.add((gene,"NOT"))


        subTerms = set()
        test = Graph()
        length = 0

        copyGraph = self.createDiGraphCopy(geneTuples)

        queue = []
        leafs = set()
        for node in self.nodes():
            if len(self.edges(node)) == 0:
                leafs.add(node)
                for parent in self.predecessors(node):
                    heappush(queue, (self.edge[parent][node]['weight'], (parent, node)))
                    
        while len(queue) > 0:
            edge = heappop(queue)

            loss = copyGraph.node[edge[1][0]]['infoLoss'] + edge[0] + copyGraph.node[edge[1][1]]['infoLoss']
            mergedNodes = copyGraph.node[edge[1][0]]['mergeCount'] + 1 + copyGraph.node[edge[1][1]]['mergeCount']
            mergedGenes = len(copyGraph.node[edge[1][0]]['mergeGene'].union(copyGraph.node[edge[1][1]]['mergeGene']).union(copyGraph.node[edge[1][1]]['gene']).union(copyGraph.node[edge[1][0]]['gene']))
            genesCount = len(copyGraph.node[edge[1][1]]['mergeGene'].union(copyGraph.node[edge[1][1]]['gene']))

            if mergedGenes <= maxMergedGeneCount and self.getLevel(edge[1][0]) >= minLevel:
                subNodes = descendantDict[edge[1][0]]
                subGraph = undirected.subgraph(subNodes.intersection(self.nodes())).copy()

                subRootGeneTuples = copyGraph.node[edge[1][0]]['mergeGene'].union(copyGraph.node[edge[1][1]]['mergeGene']).union(copyGraph.node[edge[1][1]]['gene']).union(copyGraph.node[edge[1][0]]['gene']).intersection(geneTuples)
                subRootGenes = set()
                for geneGroup in list(subRootGeneTuples):
                    subRootGenes.add(geneGroup[0])

                subTerms = set()
                for gene in list(subRootGenes):
                    subTerms.update(self.getNodesByGene(gene))
                subTerms.intersection_update(subGraph.nodes())

                subGraph = self.augmentGraph(subGraph, list(subTerms), geneTuples, fifth)
                test = make_steiner_tree(subGraph, list(subTerms))
                length = 0
                for subEdge in test.edges():
                    length += test.edge[subEdge[0]][subEdge[1]]['weight']

                prob = calcProb(length, mergedGenes, model)
            else:
                prob = 1 + maxProb

            if prob < maxProb or genesCount <= minGeneAutoMerge:
                parent = edge[1][0]
                predecessors = copyGraph.predecessors(edge[1][1])
                copyGraph.node[edge[1][0]]['mergeGene'].update(copyGraph.node[edge[1][1]]['mergeGene'])
                copyGraph.node[edge[1][0]]['mergeGene'].update(copyGraph.node[edge[1][1]]['gene'])
                copyGraph.node[edge[1][0]]['mergePMID'].update(copyGraph.node[edge[1][1]]['mergeGene'])
                copyGraph.node[edge[1][0]]['mergePMID'].update(copyGraph.node[edge[1][1]]['pmid'])
                copyGraph.remove_edge(edge[1][0], edge[1][1])
                copyGraph.remove_node(edge[1][1])
                
                for pred in predecessors:
                    if len(copyGraph.edges(pred)) == 0:
                        copyGraph = delEmptyNode(copyGraph, pred)

                copyGraph.node[parent]['infoLoss'] = length
                copyGraph.node[parent]['mergeCount'] = mergedNodes

                queue = []
                leafs = set()
                for node in copyGraph.nodes():
                    if len(copyGraph.edges(node)) == 0:
                        leafs.add(node)
                        for parent in copyGraph.predecessors(node):
                            heappush(queue, (copyGraph.edge[parent][node]['weight'], (parent, node)))

        nodes = self.nodes()

        for node in nodes:
            if node not in copyGraph:
                self.remove_node(node)
            else:
                self.node[node]['data'].addMergedGenes(copyGraph.node[node]['mergeGene'])
                self.node[node]['data'].setInfoLoss(copyGraph.node[node]['infoLoss'])
                self.node[node]['data'].setMergedCount(copyGraph.node[node]['mergeCount'])
                self.node[node]['data'].setMergedPMIDs(copyGraph.node[node]['mergePMID'])
        return self, leafs, copyGraph

    ## Creates a DiGraph copy of the graph containing only the nodes that annotate, or ancestors of nodes that annotate, the genes of interest
    # @param geneTuples The list of genes of interest
    # /return copyGraph The created DiGraph copy
    def createDiGraphCopy(self, geneTuples):
        copyGraph = DiGraph()
        for node in self:
            if len(self.getGenesByNode(node).intersection(geneTuples)) > 0 or len(self.getPropagatedGenesByNode(node).intersection(geneTuples)) > 0:
                copyGraph.add_node(node, {'gene':set(self.getGenesByNode(node).intersection(geneTuples)), 'propGene':set(self.getPropagatedGenesByNode(node).intersection(geneTuples)), 'pmid': set(self.getPubMedByNode(node)), 'mergeGene': set(), 'mergePMID': set(), 'mergeCount': 0, 'infoLoss': 0})

        for i in copyGraph:
            for edge in self.edge[i]:
                if edge in copyGraph:
                    copyGraph.add_edge(i, edge, self.edge[i][edge])

        return copyGraph

    ## Adds edges between nodes that share common genes, with the edge weight depending on the number of genes shared
    # @param subGraph The subgraph of the graph that will be augmented
    # @param termList The list of nodes within the subgraph to check for shared genes
    # @param geneTuples The genes that re of interest
    # @param fifth The edge weight that falls within the fifth percentile
    # /return subGraph The augmented subgraph
    def augmentGraph(self, subGraph, termList, geneTuples, fifth):
        newEdges = list()
        overlap = 0
        for first in range(0, len(termList)):
            for second in range(first+1, len(termList)):
                overlap = len(self.getGenesByNode(termList[first]).intersection(self.getGenesByNode(termList[second])).intersection(geneTuples))
                if overlap > 0:
                    newEdges.append([termList[first], termList[second], float(fifth)/overlap])

        for edge in newEdges:
            if edge[0] not in subGraph.edge[edge[1]]:
                subGraph.add_edge(edge[0], edge[1], {'weight':edge[2]})
            else:
                if subGraph.edge[edge[0]][edge[1]]['weight'] > edge[2]:
                    subGraph.edge[edge[0]][edge[1]]['weight'] = edge[2]
        return subGraph

    ## Checks and removes the node if it does not contain any directly associated or merged genes, also checks if any newly created leaf nodes are empty and removes them if they are
    # @param node The node that will be checked
    # /return graph The updated version of the inputted graph
    def removeEmptyNode(self, node):
        if len(self.node[node]['data'].getMergedGenes()) == 0 and len(self.getGenesByNode(node)) == 0:
            predecessors = self.predecessors(node)
            self.remove_node(node)
            for pred in predecessors:
                if len(self.edges(pred)) == 0:
                    self.removeEmptyNode(pred)


    def mergeGraphProb(self, model, maxProb=0.05, maxMergedGeneCount=200, minGeneAutoMerge=5):
        queue = []
        leafs = set()
        for node in self.nodes():
            if len(self.edges(node)) == 0:
                leafs.add(node)
                for parent in self.predecessors(node):
                    loss = self.node[parent]['data'].getInfoLoss() + self.edge[parent][node]['weight'] + self.node[node]['data'].getInfoLoss()
                    mergedGenes = len(self.node[parent]['data'].getMergedGenes().union(self.node[node]['data'].getMergedGenes()).union(self.getGenesByNode(node)).union(self.getGenesByNode(parent)))
                    prob = calcProb(loss, mergedGenes, model)
                    heappush(queue, (prob, (parent, node)))

        while len(queue) > 0:
            edge = heappop(queue)
            #Calculates information loss
            loss = self.node[edge[1][0]]['data'].getInfoLoss() + edge[0] + self.node[edge[1][1]]['data'].getInfoLoss()
            mergedNodes = self.node[edge[1][0]]['data'].getMergedCount() + 1 + self.node[edge[1][1]]['data'].getMergedCount() + 1
            mergedGenes = len(self.node[edge[1][0]]['data'].getMergedGenes().union(self.node[edge[1][1]]['data'].getMergedGenes()).union(self.getGenesByNode(edge[1][1])).union(self.getGenesByNode(edge[1][0])))
            genes = len(self.node[edge[1][1]]['data'].getMergedGenes().union(self.getGenesByNode(edge[1][1])))

            #Obtains probability
            if mergedGenes <= maxMergedGeneCount:
                prob = edge[0]
            else:
                prob = 1 + maxProb

            #Merges if probability is less than maxProb
            if prob < maxProb or genes <= minGeneAutoMerge:
                predecessors = self.predecessors(edge[1][1])
                self = mergeEdge(self, (edge[1][0], edge[1][1]))
                self.remove_node(edge[1][1])
                leafs.remove(edge[1][1])

                for pred in predecessors:
                    if len(self.edges(pred)) == 0:
                        self.removeEmptyNode(pred)

                queue = []
                for node in self.nodes():
                    if len(self.edges(node)) == 0:
                        leafs.add(node)
                        for parent in self.predecessors(node):
                            loss = self.node[parent]['data'].getInfoLoss() + self.edge[parent][node]['weight'] + self.node[node]['data'].getInfoLoss()
                            mergedGenes = len(self.node[parent]['data'].getMergedGenes().union(self.node[node]['data'].getMergedGenes()).union(self.getGenesByNode(node)))
                            if mergedGenes > maxMergedGeneCount:
                                prob = 1
                            else:
                                prob = calcProb(loss, mergedGenes, model)
                            heappush(queue, (prob, (parent, node)))
                
        return self, leafs
                
