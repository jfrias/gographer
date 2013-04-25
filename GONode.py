from utils import *

class GONode():
    def __init__ (self, goid=None, namespace=None, parents=None, obsolete=False,
                  name=None, description=None, genes=set(), propGenes=set(),
                  pmids=set(), propPmids=set(), wordVector=dict(), descendantCount=None,
                  mergedGenes=set(), mergedPmids=set(), mergedCount=0, infoLoss=0):
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
        self.descendantCount = descendantCount
        self.mergedGenes = mergedGenes
        self.mergedPmids = mergedPmids
        self.mergedCount = mergedCount
        self.infoLoss = infoLoss


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
    # @param    pmid    A list where each entry is a PubMed ID tuple, where it's the PubMed ID number followed by the qualifier
    def addPMIDs(self, pmid):
        self.pmids.update(pmid)

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
        self.propPmids.update(pmids)

    ##Set the genes associated with the node
    # @param    genes   The set containing the genes that are associated with the node.
    #                   The genes are in tuples where it's the gene ID followed by the qualifier.
    def setGenes (self, genes):
        self.genes = genes

    ##Returns the genes associated with the node
    def getGenes (self):
        return self.genes

    ##Adds list containing one or more gene tuples to the existing set of genes
    # @param    genes    A list where each entry is a gene tuple, where it's the gene ID followed by the qualifier
    def addGenes(self, genes):
        self.genes.update(genes)

    ##Set the propagated genes associated with the node or its descendants
    # @param    propGenes   The set containing the propagated PubMed IDs that are associated with the node or its descendants.
    #                   The PubMed IDs are in tuples where it's the PubMed ID number followed by the qualifier.
    def setPropagatedGenes (self, propGenes):
        self.propGenes = propGenes

    ##Returns the propagated PubMed IDs associated with the node or its descendants
    def getPropagatedGenes (self):
        return self.propGenes

    ##Adds list containing PubMed ID tuples to the existing set of propagated PubMed IDs
    # @param    pmids    A list where each entry is a PubMed ID tuple, where it's the PubMed ID number followed by the qualifier
    def addPropagatedGenes(self, genes):
        self.propGenes.update(genes)

    ##Calculates the word vector and stores this information
    # @param    corpus  The corpus that contains the information on the PubMed article
    # @param    pmids   A list of PubMed ID tuples to be added to the word vector, where it's the PubMed ID number followed by the qualifier
    #                   The propagated PubMed IDs will be used if none is given
    # @param    tokenizer   The tokenizer function that will be used on the text, a simple tokenizer is used if none is given.
    #                       Should take a string as an input, and outputs a string with words that are lower case and separated by a space
    # @param    stemmer The stemmer function that will be used to stem the words, the porter stemmer is used if none is given
    #                   Takes a word as an input and reports a stemmed word as an output.
    # @param    stopwords   A list of stop words that will not be included in the word vector, either as a list or a StopwordList.
    #                       An empty list is used if no stop word list is given.
    def calculateWordVector(self, corpus, tokenizer=Tokenizer().tokenize_word,
                            stemmer=PorterStemmer().stem, stopwords=[], pmids=None):
        wordVector = {}

        if self.getDescription():
            words = tokenizer(self.getDescription()).split(' ')
            for word in words:
                if word in stopwords:
                    continue
                word = stemmer(word, 0, len(word)-1)
                if word in wordVector:
                    wordVector[word] += 1
                else:
                    wordVector[word] = 1
        if self.getName():
            words = tokenizer(self.getName()).split(' ')
            for word in words:
                if word in stopwords:
                    continue
                word = stemmer(word, 0, len(word)-1)
                if word in wordVector:
                    wordVector[word] += 1
                else:
                    wordVector[word] = 1

        if len(self.getPropagatedPMIDs()) > 0:
            if not corpus:
                print "No corpus is given, so a word vector can not be calculated."

            else:
                if not pmids:
                    pmids = self.getPropagatedPMIDs()

                for pmid, qualifier in pmids:
                    #Checks to make sure the PubMed article is in the corpus
                    #and the qualifier does not contain the word 'NOT'
                    if pmid in corpus.docs and "NOT" not in qualifier:
                        #Adds each word to the word vector
                        pmidWordVector = corpus.docs[pmid].getWordVector(tokenizer, stemmer, stopwords)
                        for word in pmidWordVector:
                            #if word not in stopwords:
                            if word in wordVector:
                                wordVector[word] += pmidWordVector[word]
                            else:
                                wordVector[word] = pmidWordVector[word]
        self.wordVector = wordVector
    
    ##Returns the word vector for the node calculated using propagated PubMed IDs
    # @param    corpus  The corpus that contains the information on the PubMed article
    def getWordVector(self, corpus=None, stopwords=[]):
        #Calculates the word vector if the current word vector is empty
        if len(self.wordVector) == 0:
##            if len(self.getPropagatedPMIDs()) > 0:
            self.calculateWordVector(corpus, stopwords=stopwords)
        return self.wordVector

    ##Returns the number of descendants of the current node
    def getDescendantCount(self):
        return self.descendantCount

    
    ##Set the genes associated with the nodes that have been merged into the current node
    # @param    mergedGenes   The set containing the genes that are associated with the merged nodes.
    #                   The genes are in tuples where it's the gene ID followed by the qualifier.
    def setMergedGenes (self, mergedGenes):
        self.mergedGenes = mergedGenes

    ##Returns the genes associated with the nodes that have been merged into the current node
    def getMergedGenes (self):
        return self.mergedGenes

    ##Adds list containing one or more gene tuples to the existing set of merged genes
    # @param    mergedGenes    A list where each entry is a gene tuple, where it's the gene ID followed by the qualifier
    def addMergedGenes(self, mergedGenes):
        self.mergedGenes.update(mergedGenes)


    ##Set the PubMed IDs associated with the nodes that have been merged into the current node
    # @param    mergedPmids   The set containing the PubMed IDs that are associated with the merged nodes.
    #                   The PubMed IDs are in tuples where it's the PubMed ID number followed by the qualifier.
    def setMergedPMIDs (self, mergedPmids):
        self.mergedPmids = mergedPmids

    ##Returns the PubMed IDs associated with the nodes that have been merged into the current node
    def getMergedPMIDs (self):
        return self.mergedPmids

    ##Adds list containing one or more PMID tuples to the existing set of merged PMIDs
    # @param    mergedPmids    A list where each entry is a PubMed ID tuple, where it's the PubMed ID followed by the qualifier
    def addMergedPMIDs(self, mergedPmids):
        self.mergedPmids.update(mergedPmids)

    ##Set the count of the number of nodes that have been merged into the current node
    # @param    mergedCount   The integer count of the number of nodes that have been merged into the current node
    def setMergedCount (self, mergedCount):
        self.mergedCount = mergedCount

    ##Returns the PubMed IDs associated with the nodes that have been merged into the current node
    def getMergedCount (self):
        return self.mergedCount

    ##Adds list containing one or more PMID tuples to the existing set of merged PMIDs
    # @param    mergedPmids    A list where each entry is a PubMed ID tuple, where it's the PubMed ID followed by the qualifier
    def addMergedCount(self, mergedCount=1):
        self.mergedCount += mergedCount

    ##Set the PubMed IDs associated with the nodes that have been merged into the current node
    # @param    infoLoss   The set containing the PubMed IDs that are associated with the merged nodes.
    #                   The PubMed IDs are in tuples where it's the PubMed ID number followed by the qualifier.
    def setInfoLoss (self, infoLoss):
        self.infoLoss = infoLoss

    ##Returns the PubMed IDs associated with the nodes that have been merged into the current node
    def getInfoLoss (self):
        return self.infoLoss

    ##Adds list containing one or more PMID tuples to the existing set of merged PMIDs
    # @param    infoLoss    A list where each entry is a PubMed ID tuple, where it's the PubMed ID followed by the qualifier
    def addInfoLoss(self, infoLoss):
        self.infoLoss += infoLoss
        
