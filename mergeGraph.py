from calcModel import *
from heapq import *

## Merges the graph until the number of leaves left in the graph is less than or equal to the inputted leafCount
# @param graph A weighted GOGenePubmedGraph that will be merged
# @param leafCount The number of leaves that should be left in the graph when the function stops running
# /return graph The updated version fo the inputted graph 
def mergeGraph(graph, leafCount):
    queue = []
    leafs = set()
    for node in graph.nodes():
        if len(graph.edges(node)) == 0:
            leafs.add(node)
            for parent in graph.predecessors(node):
                heappush(queue, (graph.edge[parent][node]['weight'], (parent, node)))
    while len(leafs) > leafCount:
        edge = heappop(queue)
        graph = mergeEdge(graph, (edge[1][0], edge[1][1]))
        if len(graph.predecessors(edge[1][1])) == 0:
            graph.remove_node(edge[1][1])
            leafs.remove(edge[1][1])
        if len(graph.edges(edge[1][0])) == 0:
            leafs.add(edge[1][0])
            for parent in graph.predecessors(edge[1][0]):
                heappush(queue, (graph.edge[parent][edge[1][0]]['weight'], (parent, edge[1][0])))
    return graph

## Tracks the amount of information that would be lost when removing the given edge from the graph prior to removing the edge
# @param graph A weighted GOGenePubmedGraph from which to remove the edge
# @param edge The edge that is to be removed from the graph
# /return graph The updated version of the inputted graph
def mergeEdge(graph, edge):
    graph.node[edge[0]]['data'].addMergedCount(graph.node[edge[1]]['data'].getMergedCount()+1)
    graph.node[edge[0]]['data'].addMergedGenes(graph.node[edge[1]]['data'].getMergedGenes())
    graph.node[edge[0]]['data'].addMergedGenes(graph.getGenesByNode(edge[1]))
    graph.node[edge[0]]['data'].addMergedPMIDs(graph.node[edge[1]]['data'].getMergedPMIDs())
    graph.node[edge[0]]['data'].addMergedPMIDs(graph.getPubMedByNode(edge[1]))
    graph.node[edge[0]]['data'].addInfoLoss(graph.edge[edge[0]][edge[1]]['weight'])
    graph.node[edge[0]]['data'].addInfoLoss(graph.node[edge[1]]['data'].getInfoLoss())
    graph.remove_edge(edge[0], edge[1])
    return graph

def checkMerge(loss, nodeCount, model, maxProb=0.05):
    merge = False
    prob = calcProb(loss, nodeCount, model)
    if prob >= maxProb:
        merge = False
    else:
        merge = True
    return merge

## Merges the graph, and checks before each merge to ensure the probability of associated information loss is less than the set max probability
# @param graph A weighted GOGenePubmedGraph that will be merged
# @param model The model from which to calculate the probability
# @param maxProb The max allowed probaiblity for a merging to occur, the default value is 0.05
# @param maxMergedGeneCount The max allowed number of merged genes for any given node, merging stops at that node if the max is reached. The default value is 200
# /return graph The updated version fo the inputted graph
def mergeGraphCheck(graph, model, maxProb=0.05, maxMergedGeneCount=200, minGeneAutoMerge=5):
    queue = []
    leafs = set()
    for node in graph.nodes():
        if len(graph.edges(node)) == 0:
            leafs.add(node)
            for parent in graph.predecessors(node):
                heappush(queue, (graph.edge[parent][node]['weight'], (parent, node)))

    while len(queue) > 0:
        edge = heappop(queue)
        #Calculates information loss
        loss = graph.node[edge[1][0]]['data'].getInfoLoss() + edge[0] + graph.node[edge[1][1]]['data'].getInfoLoss()
        mergedNodes = graph.node[edge[1][0]]['data'].getMergedCount() + 1 + graph.node[edge[1][1]]['data'].getMergedCount() + 1
        mergedGenes = len(graph.node[edge[1][0]]['data'].getMergedGenes().union(graph.node[edge[1][1]]['data'].getMergedGenes()).union(graph.getGenesByNode(edge[1][1])))
        genes = len(graph.node[edge[1][1]]['data'].getMergedGenes().union(graph.getGenesByNode(edge[1][1])))

        #Obtains probability
        if mergedGenes <= maxMergedGeneCount:
	    prob = calcProb(loss, mergedGenes, model)
	else:
	    prob = 1 + maxProb

	#Merges if probability is less than maxProb
        if prob < maxProb or genes <= minGeneAutoMerge:
            predecessors = graph.predecessors(edge[1][1])
            graph = mergeEdge(graph, (edge[1][0], edge[1][1]))
            graph.remove_node(edge[1][1])
            leafs.remove(edge[1][1])

            for pred in predecessors:
                if len(graph.edges(pred)) == 0:
                    graph = removeEmptyNode(graph, pred)

	    queue = []
	    for node in graph.nodes():
	        if len(graph.edges(node)) == 0:
		    leafs.add(node)
		    for parent in graph.predecessors(node):
		        heappush(queue, (graph.edge[parent][node]['weight'], (parent, node)))
##            ##Add any new leaves and associated edges to the queue and leaf list
##            for pred in predecessors:
##                if len(graph.edges(pred)) == 0:
##                    leafs.add(pred)
##                    for parent in graph.predecessors(pred):
##                        heappush(queue, (graph.edge[parent][pred]['weight'], (parent, pred)))
##
##            #Check queue to remove all edges associated with removed node
##            newQueue = []
##            for i in queue:
##                if i[1][1] != edge[1][1]:
##                    heappush(newQueue, i)
##            queue = newQueue
##            
##            #Reconsiders all previously unmerged edges to leaves
##	    for i in discarded:
##		heappush(queue, i)
##	    discarded = list()
##        else:
##	    discarded.append(edge)	    
    return graph, leafs

def mergeGraphMultCheck(graph, model, maxProb=0.05, maxMergedGeneCount=200, minGeneAutoMerge=5):
    queue = []
    leafs = set()
    for node in graph.nodes():
        if len(graph.edges(node)) == 0:
            leafs.add(node)
            for parent in graph.predecessors(node):
                geneCount = len(graph.getGenesByNode(node))
                if geneCount == 0:
                    geneCount = 1
                heappush(queue, (graph.edge[parent][node]['weight']*(geneCount+len(graph.node[node]['data'].getMergedGenes())), (parent, node)))

    while len(queue) > 0:
        edge = heappop(queue)
        #Calculates information loss
        loss = graph.node[edge[1][0]]['data'].getInfoLoss() + graph.edge[edge[1][0]][edge[1][1]]['weight'] + graph.node[edge[1][1]]['data'].getInfoLoss()
        mergedNodes = graph.node[edge[1][0]]['data'].getMergedCount() + 1 + graph.node[edge[1][1]]['data'].getMergedCount() + 1
        mergedGenes = len(graph.node[edge[1][0]]['data'].getMergedGenes().union(graph.node[edge[1][1]]['data'].getMergedGenes()).union(graph.getGenesByNode(edge[1][1])))
        genes = len(graph.node[edge[1][1]]['data'].getMergedGenes().union(graph.getGenesByNode(edge[1][1])))

        #Obtains probability
        if mergedGenes <= maxMergedGeneCount:
	    prob = calcProb(loss, mergedGenes, model)
	else:
	    prob = 1 + maxProb

	#Merges if probability is less than maxProb
        if prob < maxProb or genes <= minGeneAutoMerge:
            predecessors = graph.predecessors(edge[1][1])
            graph = mergeEdge(graph, (edge[1][0], edge[1][1]))
            graph.remove_node(edge[1][1])
            leafs.remove(edge[1][1])

            for pred in predecessors:
                if len(graph.edges(pred)) == 0:
                    graph = removeEmptyNode(graph, pred)

	    queue = []
	    for node in graph.nodes():
	        if len(graph.edges(node)) == 0:
		    leafs.add(node)
		    for parent in graph.predecessors(node):
                        geneCount = len(graph.getGenesByNode(node))
                        if geneCount == 0:
                            geneCount = 1
		        heappush(queue, (graph.edge[parent][node]['weight']*(geneCount+len(graph.node[node]['data'].getMergedGenes())), (parent, node)))
	    
    return graph, leafs

## Checks and removes the node if it does not contain any directly associated or merged genes, also checks if any newly created leaf nodes are empty and removes them if they are
# @param graph A weighted GOGenePubmedGraph that will be merged
# @param node The node that will be checked
# /return graph The updated version of the inputted graph
def removeEmptyNode(graph, node):
    if len(graph.node[node]['data'].getMergedGenes()) == 0 and len(graph.getGenesByNode(node)) == 0:
        predecessors = graph.predecessors(node)
        graph.remove_node(node)
        for pred in predecessors:
            if len(graph.edges(pred)) == 0:
                removeEmptyNode(graph, pred)
    return graph


def mergeGraphProb(graph, model, maxProb=0.05, maxMergedGeneCount=200, minGeneAutoMerge=5):
    queue = []
    leafs = set()
    for node in graph.nodes():
        if len(graph.edges(node)) == 0:
            leafs.add(node)
            for parent in graph.predecessors(node):
                loss = graph.node[parent]['data'].getInfoLoss() + graph.edge[parent][node]['weight'] + graph.node[node]['data'].getInfoLoss()
                mergedGenes = len(graph.node[parent]['data'].getMergedGenes().union(graph.node[node]['data'].getMergedGenes()).union(graph.getGenesByNode(node)))
                prob = calcProb(loss, mergedGenes, model)
                heappush(queue, (prob, (parent, node)))

    while len(queue) > 0:
        edge = heappop(queue)
        #Calculates information loss
        loss = graph.node[edge[1][0]]['data'].getInfoLoss() + edge[0] + graph.node[edge[1][1]]['data'].getInfoLoss()
        mergedNodes = graph.node[edge[1][0]]['data'].getMergedCount() + 1 + graph.node[edge[1][1]]['data'].getMergedCount() + 1
        mergedGenes = len(graph.node[edge[1][0]]['data'].getMergedGenes().union(graph.node[edge[1][1]]['data'].getMergedGenes()).union(graph.getGenesByNode(edge[1][1])))
        genes = len(graph.node[edge[1][1]]['data'].getMergedGenes().union(graph.getGenesByNode(edge[1][1])))

        #Obtains probability
        if mergedGenes > maxMergedGeneCount:
	    prob = edge[0]
	else:
	    prob = 1 + maxProb

	#Merges if probability is less than maxProb
        if prob < maxProb or genes <= minGeneAutoMerge:
            predecessors = graph.predecessors(edge[1][1])
            graph = mergeEdge(graph, (edge[1][0], edge[1][1]))
            graph.remove_node(edge[1][1])
            leafs.remove(edge[1][1])

            for pred in predecessors:
                if len(graph.edges(pred)) == 0:
                    graph = removeEmptyNode(graph, pred)

	    queue = []
	    for node in graph.nodes():
	        if len(graph.edges(node)) == 0:
		    leafs.add(node)
		    for parent in graph.predecessors(node):
                        loss = graph.node[parent]['data'].getInfoLoss() + graph.edge[parent][node]['weight'] + graph.node[node]['data'].getInfoLoss()
                        mergedGenes = len(graph.node[parent]['data'].getMergedGenes().union(graph.node[node]['data'].getMergedGenes()).union(graph.getGenesByNode(node)))
                        prob = calcProb(loss, mergedGenes, model)
		        heappush(queue, (prob, (parent, node)))
	    
    return graph, leafs
