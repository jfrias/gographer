## A Python object representation of a PubMed record.
class Document:

    ## Constructor
    def __init__(self, pmid = None, abstract = None, title = None):
        self.pmid = pmid
        self.abstract = abstract
        self.title = title
        self.wordVec = {}
        self.meshTerms = {}


    ## Add word to wordVec
    # @param word  The index of word
    # @param count The count to be added into the wordVec
    def addWord(self, word, count = 1):
        if word in self.wordVec:
            self.wordVec[word] += count
        else:
            self.wordVec[word] = count


    ## Add MeSH Term to Document
    # @param mesh The index of the MeSH term
    # @param count The count to be added into the MeSHTerm
    def addMeSH(self, mesh, count = 1):
        if mesh in self.meshTerms:
            self.meshTerms[mesh] += count
        else:
            self.meshTerms[mesh] = count


    ## Determines whether self has a MeSH Term contained within
    def hasMeSH(self):
        return len(self.meshTerms) > 0


    ## Evaluating equality
    def __eq__(self,other):
        if not isinstance(other, Document):
            raise TypeError, "other must be an instance of Document"
        if self.pmid is None or other.pmid is None:
            raise RuntimeError, "cannot compare documents when one of the PMID's is None"
        return self.pmid == other.pmid and self.abstract == other.abstract

        
    def __str__(self):
        return "<Document: %s>" % self.pmid
    

    def __repr__(self):
        return str(self)

    ## Returns the word vector of the Documented calculated on the fly
    def getWordVector(self, tokenizer, stemmer):
        wordVector = {}
        
        words = tokenizer(self.title).split(' ')
        for word in words:
            word = stemmer(word, 0, len(word)-1)
            if word in wordVector:
                wordVector[word] += 1
            else:
                wordVector[word] = 1

        words = tokenizer(self.abstract).split(' ')
        for word in words:
            word = stemmer(word, 0, len(word)-1)
            if word in wordVector:
                wordVector[word] += 1
            else:
                wordVector[word] = 1
        return wordVector
        
