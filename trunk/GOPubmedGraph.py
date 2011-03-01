## @package GOPubmedGraph

from GOGraph import GOGraph

class GOPubmedGraph(GOGraph):
    ## Create a pubmed graph from a GOGraph
    # @param    gograph A GOGraph to base this graph off of
    # @param    assoc   The file containing gene and pubmed association information
    # @param    excludeEvidence A list of the evidence codes that should be ignored
    def __init__(self, gograph, assoc=None, excludeEvidence =[]):
        if not gograph.__class__.__name__ is 'GOGraph':
            msg =  "You did not give an instance of GOGraph to GOProteinGraph (it's a %s):" % gograph.__class__.__name__ \
                  + " any weights or directions will be ignored."
            
        GOGraph.__init__(self)
        self.add_edges_from(gograph.edges_iter(data=True))
        self.add_nodes_from(gograph.nodes_iter(data=True))
        self.pubmedToNode = dict()
        self.excludeEvidence = excludeEvidence

        #Adds a 'pubmd' field for all nodes in the graph
        for node in self.nodes():
            self.node[node]['pubmed'] = list()

        if assoc != None:
            self.parseAssocFile(assoc)

    ## Parses the given association file and adds the pubmed information to the appropriate nodes
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

                    #Checks to make sure the evidence type isn't one that was listed to be excluded
                    if fields[6] in self.excludeEvidence:
                        continue

                    # Pmid's that start w/ P may bes provincial
                    references = filter(lambda x: x.startswith("PMID:") and not x.startswith("PMID:P"), fields[5].split('|'))

                    if(len(references) > 0):
                        for ref in references:
                            #Stores both the id and the qualifier in the node's pubmed list
                            self.node[fields[4]]['pubmed'].append((ref[5:],fields[3]))

                            #Adds the pubmed to node association to the dictionary only if it isn't already in dictionary
                            if ref[5:] not in self.pubmedToNode:
                                self.pubmedToNode[ref[5:]] = [fields[4]]
                            else:
                                if fields[4] not in self.pubmedToNode[ref[5:]]:
                                    self.pubmedToNode[ref[5:]].append(fields[4])
            f.close
        except:
            print "Could not parse association file %s" % (assoc)

    ## Applies a weight to the graph
    def weight(self):
        raise NotImplementedError

    ## Returns a GOGenePubmedGraph version of itself
    # @param    gogenegraph The GOGeneGraph that this graph is to be merged with
    def toGOGenePubmedGraph(self):
        raise NotImplementedError

    ## Returns the associated PubMed IDs for a given node
    # @param    goid    The GOID of the node to retrieve the associated PubMed IDs from
    def getPubMedByNode(self, goid):
        if goid in self:
            return self.node[goid]['pubmed']
        else:
            print "Given GOID is not in the graph"
            raise

    ## Returns the associated nodes for a given PubMed ID
    # @param    pubmed  The PubMed ID for which to retrieve the associated nodes:
    def getNodesByPubMed(self, pubmed):
        if pubmed in self.pubmedToNode:
            return self.pubmedToNode[pubmed]
        else:
            print "Given pubmed id is not in the graph"
            raise

