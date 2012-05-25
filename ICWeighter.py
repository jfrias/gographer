from WeightingInterface import WeightingInterface
from networkx import topological_sort
import math

class ICWeighter(WeightingInterface):
    def __init__(self, graph=None, smoother=0.001):
        self.graphID = None
        self.totalAnnotation = 0
        self.smoother = smoother
        if graph:
            self.storeTotalAnnotation(graph)
    
    def makeWeighted(self, child, parent, graph):
        if id(graph) != self.graphID:
            self.storeTotalAnnotation(graph)
        return self.calcWeightIC(child, parent, graph)
            
    def storeTotalAnnotation(self, graph):
        self.graphID = id(graph)
        totalAnnotation = 0
        for node in graph:
            totalAnnotation += len(graph.getPropagatedPubMedByNode(node).union(graph.getPubMedByNode(node)))
        self.totalAnnotation = totalAnnotation

    def setSmoother(self, smoother):
        self.smoother = smoother

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
        
        
