from WeightingInterface import WeightingInterface
from IBWeighter import *
import math

class GeneIBWeighter(WeightingInterface):

	## Add weights to the connections from parents to self to children
	# @param	child The child node
	# @param	parent The parent node
	# @param	graph The entire graph on which these nodes exist
    def makeWeighted(self, child, parent, graph):
        return self.lossIB(child, parent, graph) + self.lossProtein(child, parent, graph)

	## IB loss calculation
	# @param	child The child node
	# @param	parent The parent node
	# @param	graph The entire graph on which these nodes exist
    def lossIB(self, child, parent, graph):
        childNode = graph.node[child]['data']
        parentNode = graph.node[parent]['data']

        
        tChild = float(graph.getDescendantCount(child)+1)
        tParent = float(graph.getDescendantCount(parent)+1)
        tRoot = float(graph.number_of_nodes())

        pti = tChild/tRoot
        pii = tChild/tParent
        distanceKL = self.calcWeightKL(childNode, parentNode)

        return pti*pii*distanceKL
	## Calculate the KL Weight
	# @param	n1
	# @param	n2
	# @param	smoother Always set to 0.001 automatically
    def calcWeightKL(self, n1, n2, smoother=0.001):
        v1 = n1.getWordVector().copy()
        v2 = n2.getWordVector().copy()

        total1 = 0.0
        total2 = 0.0

        for w in v1:
            total1 += v1[w]
        for w in v2:
            total2 += v2[w]

        for w in v1:
            if w not in v2:
                v2[w] = smoother
        for w in v2:
            if w not in v1:
                v1[w] = smoother

        distance = 0
        for w in v1:
            distance += v1[w]/total1 * math.log((v1[w]/total1)/(v2[w]/total2))
        return distance

	## Protein loss calculation
	# @param	child The child node
	# @param	parent The parent node
	# @param	graph The entire graph on which these nodes exist
    def lossProtein(self, child, parent, graph):
        childProteins = float(len(graph.node[child]['data'].getPropagatedGenes()))
        parentProteins = float(len(graph.node[parent]['data'].getPropagatedGenes()))
        totalProteins = float(len(graph.geneToNode))
        return math.log(parentProteins/totalProteins) - math.log(childProteins/totalProteins)
