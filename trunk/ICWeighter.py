from WeightingInterface import WeightingInterface
from networkx import topological_sort
import math

class ICWeighter(WeightingInterface):
    ## Constructor
	# @param	graph Always set to None
	# @param	smoother Always set to 0.001
    def __init__(self, graph=None, smoother=0.001):
        self.graphID = None
        self.totalAnnotation = 0
        self.smoother = smoother
        if graph:
            self.storeTotalAnnotation(graph)
    
	## Add weights to a section of the graph
	# @param	child The child node
	# @param	parent The parent node
	# @param	graph The entire graph which contains the child and parent nodes
    def makeWeighted(self, child, parent, graph):
        if id(graph) != self.graphID:
            self.storeTotalAnnotation(graph)
        return self.calcWeightIC(child, parent, graph)
    
	## Find and save the total annotation
	# @param	graph The entire graph
    def storeTotalAnnotation(self, graph):
        self.graphID = id(graph)
        totalAnnotation = 0
        for node in graph:
            totalAnnotation += len(graph.getPropagatedPubMedByNode(node).union(graph.getPubMedByNode(node)))
        self.totalAnnotation = totalAnnotation

	## Sets the value of the smoother
	# @param	smoother Value which smoother will be set to 
    def setSmoother(self, smoother):
        self.smoother = smoother

	## IC Weight calculater
	# @param	child The child node
	# @param	parent The parent node
	# @param	graph The entire graph which contains the child and parent nodes
    def calcWeightIC(self, child, parent, graph):
        childAnnotation = len(graph.getPropagatedPubMedByNode(child).union(graph.getPubMedByNode(child)))
        parentAnnotation = len(graph.getPropagatedPubMedByNode(parent).union(graph.getPubMedByNode(parent)))
        if childAnnotation == 0:
            childAnnotation = self.smoother
        if parentAnnotation == 0:
            parentAnnotation = self.smoother
        childIC = -math.log(float(childAnnotation)/self.totalAnnotation)
        parentIC = -math.log(float(parentAnnotation)/self.totalAnnotation)
        return abs(parentIC-childIC)
        
        
