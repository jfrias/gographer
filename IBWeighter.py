from WeightingInterface import WeightingInterface
from networkx import topological_sort
import math

class IBWeighter(WeightingInterface):

    def makeWeighted(self, child, parent, graph):
        return self.lossIB(child, parent, graph)

    ##Calculates the semantic distance between a child and a parent node using the information bottleneck method
    # @param    child   The GO ID of the child node
    # @param    parent  The GO ID of the parent node
    # @param    graph   The graph that contains the nodes
    def lossIB(self, child, parent, graph):
        childNode = graph.node[child]['data']
        parentNode = graph.node[parent]['data']

        
        tChild = float(graph.getDescendantCount(child)+1)
        tParent = float(graph.getDescendantCount(parent)+1)
        tRoot = float(graph.getDescendantCount(topological_sort(graph)[0])+1)

        pti = tChild/tRoot
        pii = tChild/tParent
        distanceKL = self.calcWeightKL(childNode, parentNode)

        return pti*pii*distanceKL

    ##Calculates the Kullback-Leibler weight for the two given nodes
    # @param    n1  The first node to be used in the calculation
    # @param    n2  The second node to be used in the calculation
    # @param    smoother    The smoothing value to be used if any words are missing in either word vector
    def calcWeightKL(self, n1, n2, smoother=0.001):
        v1 = n1.getWordVector().copy()
        v2 = n2.getWordVector().copy()

        total1 = 0.0
        total2 = 0.0

        for w in v1:
            if w not in v2:
                v2[w] = smoother
        for w in v2:
            if w not in v1:
                v1[w] = smoother

        for w in v1:
            total1 += v1[w]
        for w in v2:
            total2 += v2[w]

        distance = 0
        for w in v1:
            distance += v1[w]/total1 * math.log((v1[w]/total1)/(v2[w]/total2))
        return distance
