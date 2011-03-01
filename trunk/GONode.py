class GONode():
    def __init__ (self, goid=None, namespace=None, parents=None, obsolete=False, name=None, description=None):
        self.goid = goid
        self.namespace = namespace
        self.parents = parents
        self.obsolete = obsolete
        self.name = name
        self.description = description

    def setGOID (self, goid):
        self.goid = goid

    def setNamespace (self, namespace):
        self.namespace = namespace

    def setParents (self, parents):
        self.parents = parents

    def setObsolete (self, obsolete):
        self.obsolete = obsolete

    def setName (self, name):
        self.name = name

    def setDescription (self, description):
        self.description = description