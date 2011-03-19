class GONode():
    def __init__ (self, goid=None, namespace=None, parents=None, obsolete=False,
                  name=None, description=None, genes=set(), propGenes=set(),
                  pmids=set(), propPmids=set(), wordVector=dict()):
        self.goid = goid
        self.namespace = namespace
        self.parents = parents
        self.obsolete = obsolete
        self.name = name
        self.description = description
        self.genes = genes
        self.propGenes = propGenes
        self.pmids = pmids
        self.propPmids = propPmids
        self.wordVector = wordVector

    ##Set the GO ID of the node
    # @param    goid    The GO ID that will be assigned to the node
    def setGOID (self, goid):
        self.goid = goid

    ##Returns the GO ID of the node
    def getGOID (self):
        return self.goid

    ##Set the namespace of the node
    # @param    namespace   The namespace that will be assigned to the node 
    def setNamespace (self, namespace):
        self.namespace = namespace

    ##Returns the namespace of the node
    def getNamespace (self):
        return self.namespace

    ##Set the parents of the node
    # @param    parents  The parents that will be assigned to the node
    def setParents (self, parents):
        self.parents = parents
        
    ##Returns the parents of the node
    def getParents (self):
        return self.parents

    ##Set the obsolete status of the node
    # @param    obsolete    The obsolete status that will be assigned to the node
    def setObsolete (self, obsolete):
        self.obsolete = obsolete

    ##Returns the obsolete status of the node (whether or not the term is obsolete)
    def getObsolete (self):
        return self.obsolete

    ##Set the name of the node
    # @param    name    The name that will be assigned to the node
    def setName (self, name):
        self.name = name

    ##Returns the name of the node
    def getName (self):
        return self.name
    
    ##Set the description of the node
    # @param    description The description that will be assigned to the node
    def setDescription (self, description):
        self.description = description

    ##Returns the description of the node
    def getDescription (self):
        return self.description

    ##Set the PubMed IDs associated with the node
    # @param    pmids   The set containing the PubMed IDs that are associated with the node.
    #                   The PubMed IDs are in tuples where it's the PubMed ID number followed by the qualifier.
    def setPMIDs (self, pmids):
        self.pmids = pmids

    ##Returns the PubMed IDs associated with the node
    def getPMIDs (self):
        return self.pmids

    ##Adds list containing one or more PubMed ID tuples to the existing set of PubMed IDs
    # @param    pmid    A list where each entry is a PubMed IDs tuple, where it's the PubMed ID number followed by the qualifier
    def addPMIDs(self, pmid):
        self.pmids = self.pmids.union(pmid)

    ##Set the propagated PubMed IDs associated with the node or its descendants
    # @param    propPmids   The set containing the propagated PubMed IDs that are associated with the node or its descendants.
    #                   The PubMed IDs are in tuples where it's the PubMed ID number followed by the qualifier.
    def setPropagatedPMIDs (self, propPmids):
        self.propPmids = propPmids

    ##Returns the propagated PubMed IDs associated with the node or its descendants
    def getPropagatedPMIDs (self):
        return self.propPmids

    ##Adds list containing PubMed ID tuples to the existing set of propagated PubMed IDs
    # @param    pmids    A list where each entry is a PubMed ID tuple, where it's the PubMed ID number followed by the qualifier
    def addPropagatedPMIDs(self, pmids):
        self.propPmids = self.propPmids.union(pmids)
        
