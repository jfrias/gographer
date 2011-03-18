class GONode():
    def __init__ (self, goid=None, namespace=None, parents=None, obsolete=False, name=None, description=None):
        self.goid = goid
        self.namespace = namespace
        self.parents = parents
        self.obsolete = obsolete
        self.name = name
        self.description = description

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
