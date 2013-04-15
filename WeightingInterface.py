class WeightingInterface:
    
	## Weight a graph
	# @param	child A node
	# @param	parent A node
	# @param	graph An unweighted GOPubmedGraph
    def makeWeighted(self, child, parent, graph):
        raise NotImplementedError
