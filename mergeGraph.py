from calcModel import *
from heapq import *

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

def mergeEdge(graph, edge):
    graph.node[edge[0]]['data'].addMergedCount(graph.node[edge[1]]['data'].getMergedCount()+1)
    graph.node[edge[0]]['data'].addMergedGenes(graph.node[edge[1]]['data'].getPropagatedGenes())
    graph.node[edge[0]]['data'].addMergedPMIDs(graph.node[edge[1]]['data'].getPropagatedPMIDs())
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

def mergeGraphCheck(graph, model, maxProb=0.05):
    queue = []
    leafs = set()
    for node in graph.nodes():
        if len(graph.edges(node)) == 0:
            leafs.add(node)
            for parent in graph.predecessors(node):
                heappush(queue, (graph.edge[parent][node]['weight'], (parent, node)))
    discarded = list()
    while len(queue) > 0:
        edge = heappop(queue)
        loss = graph.node[edge[1][0]]['data'].getInfoLoss() + edge[0] + graph.node[edge[1][1]]['data'].getInfoLoss()
        mergedNodes = graph.node[edge[1][0]]['data'].getMergedCount() + 1 + graph.node[edge[1][1]]['data'].getMergedCount() + 1
        mergedGenes = len(graph.node[edge[1][0]]['data'].getMergedGenes().union(graph.node[edge[1][1]]['data'].getMergedGenes()).union(graph.node[edge[1][1]]['data'].getPropagatedGenes()))
        if mergedGenes <= 200:
	    prob = calcProb(loss, mergedGenes, model)
	else:
	    prob = 1 + maxProb
        if prob < maxProb:
            graph = mergeEdge(graph, (edge[1][0], edge[1][1]))
            if len(graph.predecessors(edge[1][1])) == 0:
                graph.remove_node(edge[1][1])
                leafs.remove(edge[1][1])
            if len(graph.edges(edge[1][0])) == 0:
                leafs.add(edge[1][0])
                for parent in graph.predecessors(edge[1][0]):
                    heappush(queue, (graph.edge[parent][edge[1][0]]['weight'], (parent, edge[1][0])))
	    for i in discarded:
		heappush(queue, i)
	    discarded = list()
        else:
	    discarded.append(edge)
    return graph, leafs
