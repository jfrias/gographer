from xml.sax import make_parser
from parsers import *
from Tokenizer import Tokenizer
from PorterStemmer import PorterStemmer

## A collection of Document objects.  Provide functionalities for storing and retrieval of Documents.
class Corpus:

    ## Constructor
    # @param docs An initial dictionary of the form [<pmid>] = <Document>
    def __init__(self, docs=None, tokenizer=None, stemmer=None):
        self.docs = docs or {}
        self.tokenizer = tokenizer or Tokenizer().tokenize_word
        self.stemmer = stemmer or PorterStemmer().stem


    ## Initialize corpus from a PubmedArticleSet xml file
    # @param filename The location of the xml file.
    # \return A Corpus instance
    @classmethod
    def fromPubmedArticleSetFile(klass, filename):
        handler = PubmedArticleSet.parse(filename)
        return Corpus(docs = handler.docs)


    @classmethod
    def fromGOADocumentSetFile(klass, filename):
        handler = GOADocumentSet.parse(filename)
        return Corpus(docs = handler.docs)        


    ## Get all pmids in the corpus
    # \return A list of pmids in this corpus
    def pmids(self):
        return self.docs.keys()


    ## Iterator for this corpus - iterates through Document objects
    def __iter__(self):
        for doc in self.docs.values():
            yield doc


    ## Append a Document to this corpus
    # @param doc A Document to append
    def append(self, doc):
        if doc.pmid == None:
            raise Exception, "Attempt to append Document to a Corpus with a null pmid attribute."
        self.docs[doc.pmid] = doc
            


    ## Get a Document with a given PMID
    # @param pmid The PMID of the Document to get
    # \return A Document with given PMID
    def __getitem__(self, pmid):
        return self.docs[pmid]


    ## Get number of Documents in this Corpus
    # \return Number of Documents in this Corpus
    def __len__(self):
        return len(self.docs)

    ## Calculates the word vector for all Documents in this Corpus using provided tokenizer and stemmer
    # @param tokenizer The tokenizer function that will be used on the text
    # @param stemmer The stemmer function that will be used on the text
    def calculateWordVectors(self, tokenizer=None, stemmer=None):
        if tokenizer:
            self.tokenizer = tokenizer
        if stemmer:
            self.stemmer = stemmer

        for doc in iter(self):
            title = doc.title
            title = self.tokenizer(title)
            words = title.split(' ')
            for word in words:
                doc.addWord(self.stemmer(word, 0, len(word)-1))

            abstract = doc.abstract
            abstract = self.tokenizer(abstract)
            words = abstract.split(' ')
            for word in words:
                doc.addWord(self.stemmer(word, 0, len(word)-1))
        
