#!/usr/bin/env python

import urllib, sys, time, commands, os
import string, re

from mergeGraph import *
from parsers import PubmedArticleSet

## Check to make sure all the pmids given can be found in the XML file given.
# @param ids The ids of the documents to look for.
# @param filename The PubmedArticleSet XML file to look at
# \return True if all ids found (and exactly those ids), False otherwise
def verifyDocuments(ids, filename):
    handler = PubmedArticleSet.parse(filename)
    pmids = handler.docs.keys()
    
    notfound = set()
    extrafound = set()
    for pmid in pmids:
        if not pmid in ids:
            extrafound.add(pmid)
    for gid in ids:
        if not gid in pmids:
            notfound.add(gid)

    result = True
    if len(extrafound) > 0:
        #raise RuntimeWarning, "There were %i extra pmids downloaded: %s" % (len(extrafound), ",".join(extrafound))
        print "There were %i extra pmids downloaded: %s" % (len(extrafound), ",".join(extrafound))
        result = False
    if len(notfound) > 0:
        #raise RuntimeWarning, "There were %i pmids not downloaded: %s" % (len(notfound), ",".join(notfound))                
        print "There were %i pmids not downloaded: %s" % (len(notfound), ",".join(notfound))                
        result = False
    return result


## Fetch the pubmed documents given by their ids and store them in
# xml format in the given directory
# @param ids The ids to fetch
# @param directory The directory to store the files in
# @param blocksize The number to download at one time, default of 500
# @param failIfProblem End program with error if there is a problem downloading
# pubmed documents.  If False (default) only a warning message will be displayed.
def efetchByBlock(ids, directory, blocksize=None, failIfProblem=False):
    blocksize = blocksize or 500
    url = "http://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    ids = map(lambda x: str(int(x)), ids)

    size = int(float(len(ids)) / float(blocksize))
    try:
        index = 0
        while index < len(ids):
            end = min(index+blocksize, len(ids))
            filename = os.path.join(directory, "block_%i.xml" % (index / blocksize))
            if not os.path.exists(filename):                
                outfile = open(filename, 'w')
                idblock = ",".join(ids[index:end])
                mysock = urllib.urlopen("%s?db=pubmed&id=%s&retmode=xml" % (url, idblock))
                line = mysock.readline()
                while line:
                    outfile.write(line)
                    line = mysock.readline()
                mysock.close()
                outfile.close()
                time.sleep(3)

            if not verifyDocuments(ids[index:end], filename) and failIfProblem:
                print "PMID check verification failed - some documents not downloaded"
                #raise RuntimeError, "PMID check verification failed - some documents not downloaded"
            index += blocksize
    except Exception, e:
        raise RuntimeError, "Problem fetching all PMIDs."

class Tokenizer:
    
    def __init__(self):
        self.pattern_punctuation = re.compile('(\!|\.|\+|\?|\;|\,|\:|\&|\_|\~|\:|\-|\;|\$|\@|\#|\%|\*|\^|\`|\|)')
        self.pattern_numbers = re.compile('[0-9]')
        self.pattern_plurals = re.compile('\'s$')
        self.pattern_possessives = re.compile('s\'$')
        self.pattern_blanks = re.compile('[\t, ,\n,\r,\f,\v,]+')
        self.pattern_others = re.compile('(\'|\"|\/|\\\|\<|\>|\{|\}|\[|\]|\(|\))')
        
        
        self.pattern_words_and_blanks = re.compile('[\W\t\n\r\f\v| ]|\|')
        #self.pattern.word = re.compile('[a-z]+')[^a-zA-Z]

        self.leftchar = {}
        self.leftchar[' ']=0
        for i in string.ascii_letters:
            self.leftchar[i]=0
        #for i in range(10):
            #self.leftchar[str(i)]=0
        #print len(self.leftchar)

            
    def keepNumChr(self, words):
        newWords = ''
        for c in words:
            if self.leftchar.has_key(c):
                newWords +=c
        return newWords
    
                                
    def tokenize_word(self,words):
        before = words
        words = words.lower()
        words = self.pattern_plurals.sub('', words)
        words = self.pattern_possessives.sub('', words)
        words = self.keepNumChr(words)
        words = self.pattern_blanks.sub(' ',words)
        
        #words = string.strip(words)    
        return words

        
    def test(self,words):
        for c in words:    
            if not self.leftchar.has_key(c):
                self.leftchar[c] = 1

                
    def returnUnexpChar(self):
        reV = []
        for k in self.leftchar.keys():
            if self.leftchar[k] != 0:
                reV +=[k]
                
        return reV

class PorterStemmer:

    def __init__(self):
        """The main part of the stemming algorithm starts here.
        b is a buffer holding a word to be stemmed. The letters are in b[k0],
        b[k0+1] ... ending at b[k]. In fact k0 = 0 in this demo program. k is
        readjusted downwards as the stemming progresses. Zero termination is
        not in fact used in the algorithm.

        Note that only lower case sequences are stemmed. Forcing to lower case
        should be done before stem(...) is called.
        """

        self.b = ""  # buffer for word to be stemmed
        self.k = 0
        self.k0 = 0
        self.j = 0   # j is a general offset into the string

    def cons(self, i):
        """cons(i) is TRUE <=> b[i] is a consonant."""
        if self.b[i] == 'a' or self.b[i] == 'e' or self.b[i] == 'i' or self.b[i] == 'o' or self.b[i] == 'u':
            return 0
        if self.b[i] == 'y':
            if i == self.k0:
                return 1
            else:
                return (not self.cons(i - 1))
        return 1

    def m(self):
        """m() measures the number of consonant sequences between k0 and j.
        if c is a consonant sequence and v a vowel sequence, and <..>
        indicates arbitrary presence,

           <c><v>       gives 0
           <c>vc<v>     gives 1
           <c>vcvc<v>   gives 2
           <c>vcvcvc<v> gives 3
           ....
        """
        n = 0
        i = self.k0
        while 1:
            if i > self.j:
                return n
            if not self.cons(i):
                break
            i = i + 1
        i = i + 1
        while 1:
            while 1:
                if i > self.j:
                    return n
                if self.cons(i):
                    break
                i = i + 1
            i = i + 1
            n = n + 1
            while 1:
                if i > self.j:
                    return n
                if not self.cons(i):
                    break
                i = i + 1
            i = i + 1

    def vowelinstem(self):
        """vowelinstem() is TRUE <=> k0,...j contains a vowel"""
        for i in range(self.k0, self.j + 1):
            if not self.cons(i):
                return 1
        return 0

    def doublec(self, j):
        """doublec(j) is TRUE <=> j,(j-1) contain a double consonant."""
        if j < (self.k0 + 1):
            return 0
        if (self.b[j] != self.b[j-1]):
            return 0
        return self.cons(j)

    def cvc(self, i):
        """cvc(i) is TRUE <=> i-2,i-1,i has the form consonant - vowel - consonant
        and also if the second c is not w,x or y. this is used when trying to
        restore an e at the end of a short  e.g.

           cav(e), lov(e), hop(e), crim(e), but
           snow, box, tray.
        """
        if i < (self.k0 + 2) or not self.cons(i) or self.cons(i-1) or not self.cons(i-2):
            return 0
        ch = self.b[i]
        if ch == 'w' or ch == 'x' or ch == 'y':
            return 0
        return 1

    def ends(self, s):
        """ends(s) is TRUE <=> k0,...k ends with the string s."""
        length = len(s)
        if s[length - 1] != self.b[self.k]: # tiny speed-up
            return 0
        if length > (self.k - self.k0 + 1):
            return 0
        if self.b[self.k-length+1:self.k+1] != s:
            return 0
        self.j = self.k - length
        return 1

    def setto(self, s):
        """setto(s) sets (j+1),...k to the characters in the string s, readjusting k."""
        length = len(s)
        self.b = self.b[:self.j+1] + s + self.b[self.j+length+1:]
        self.k = self.j + length

    def r(self, s):
        """r(s) is used further down."""
        if self.m() > 0:
            self.setto(s)

    def step1ab(self):
        """step1ab() gets rid of plurals and -ed or -ing. e.g.

           caresses  ->  caress
           ponies    ->  poni
           ties      ->  ti
           caress    ->  caress
           cats      ->  cat

           feed      ->  feed
           agreed    ->  agree
           disabled  ->  disable

           matting   ->  mat
           mating    ->  mate
           meeting   ->  meet
           milling   ->  mill
           messing   ->  mess

           meetings  ->  meet
        """
        if self.b[self.k] == 's':
            if self.ends("sses"):
                self.k = self.k - 2
            elif self.ends("ies"):
                self.setto("i")
            elif self.b[self.k - 1] != 's':
                self.k = self.k - 1
        if self.ends("eed"):
            if self.m() > 0:
                self.k = self.k - 1
        elif (self.ends("ed") or self.ends("ing")) and self.vowelinstem():
            self.k = self.j
            if self.ends("at"):   self.setto("ate")
            elif self.ends("bl"): self.setto("ble")
            elif self.ends("iz"): self.setto("ize")
            elif self.doublec(self.k):
                self.k = self.k - 1
                ch = self.b[self.k]
                if ch == 'l' or ch == 's' or ch == 'z':
                    self.k = self.k + 1
            elif (self.m() == 1 and self.cvc(self.k)):
                self.setto("e")

    def step1c(self):
        """step1c() turns terminal y to i when there is another vowel in the stem."""
        if (self.ends("y") and self.vowelinstem()):
            self.b = self.b[:self.k] + 'i' + self.b[self.k+1:]

    def step2(self):
        """step2() maps double suffices to single ones.
        so -ization ( = -ize plus -ation) maps to -ize etc. note that the
        string before the suffix must give m() > 0.
        """
        if self.b[self.k - 1] == 'a':
            if self.ends("ational"):   self.r("ate")
            elif self.ends("tional"):  self.r("tion")
        elif self.b[self.k - 1] == 'c':
            if self.ends("enci"):      self.r("ence")
            elif self.ends("anci"):    self.r("ance")
        elif self.b[self.k - 1] == 'e':
            if self.ends("izer"):      self.r("ize")
        elif self.b[self.k - 1] == 'l':
            if self.ends("bli"):       self.r("ble") # --DEPARTURE--
            # To match the published algorithm, replace this phrase with
            #   if self.ends("abli"):      self.r("able")
            elif self.ends("alli"):    self.r("al")
            elif self.ends("entli"):   self.r("ent")
            elif self.ends("eli"):     self.r("e")
            elif self.ends("ousli"):   self.r("ous")
        elif self.b[self.k - 1] == 'o':
            if self.ends("ization"):   self.r("ize")
            elif self.ends("ation"):   self.r("ate")
            elif self.ends("ator"):    self.r("ate")
        elif self.b[self.k - 1] == 's':
            if self.ends("alism"):     self.r("al")
            elif self.ends("iveness"): self.r("ive")
            elif self.ends("fulness"): self.r("ful")
            elif self.ends("ousness"): self.r("ous")
        elif self.b[self.k - 1] == 't':
            if self.ends("aliti"):     self.r("al")
            elif self.ends("iviti"):   self.r("ive")
            elif self.ends("biliti"):  self.r("ble")
        elif self.b[self.k - 1] == 'g': # --DEPARTURE--
            if self.ends("logi"):      self.r("log")
        # To match the published algorithm, delete this phrase

    def step3(self):
        """step3() dels with -ic-, -full, -ness etc. similar strategy to step2."""
        if self.b[self.k] == 'e':
            if self.ends("icate"):     self.r("ic")
            elif self.ends("ative"):   self.r("")
            elif self.ends("alize"):   self.r("al")
        elif self.b[self.k] == 'i':
            if self.ends("iciti"):     self.r("ic")
        elif self.b[self.k] == 'l':
            if self.ends("ical"):      self.r("ic")
            elif self.ends("ful"):     self.r("")
        elif self.b[self.k] == 's':
            if self.ends("ness"):      self.r("")

    def step4(self):
        """step4() takes off -ant, -ence etc., in context <c>vcvc<v>."""
        if self.b[self.k - 1] == 'a':
            if self.ends("al"): pass
            else: return
        elif self.b[self.k - 1] == 'c':
            if self.ends("ance"): pass
            elif self.ends("ence"): pass
            else: return
        elif self.b[self.k - 1] == 'e':
            if self.ends("er"): pass
            else: return
        elif self.b[self.k - 1] == 'i':
            if self.ends("ic"): pass
            else: return
        elif self.b[self.k - 1] == 'l':
            if self.ends("able"): pass
            elif self.ends("ible"): pass
            else: return
        elif self.b[self.k - 1] == 'n':
            if self.ends("ant"): pass
            elif self.ends("ement"): pass
            elif self.ends("ment"): pass
            elif self.ends("ent"): pass
            else: return
        elif self.b[self.k - 1] == 'o':
            if self.ends("ion") and (self.b[self.j] == 's' or self.b[self.j] == 't'): pass
            elif self.ends("ou"): pass
            # takes care of -ous
            else: return
        elif self.b[self.k - 1] == 's':
            if self.ends("ism"): pass
            else: return
        elif self.b[self.k - 1] == 't':
            if self.ends("ate"): pass
            elif self.ends("iti"): pass
            else: return
        elif self.b[self.k - 1] == 'u':
            if self.ends("ous"): pass
            else: return
        elif self.b[self.k - 1] == 'v':
            if self.ends("ive"): pass
            else: return
        elif self.b[self.k - 1] == 'z':
            if self.ends("ize"): pass
            else: return
        else:
            return
        if self.m() > 1:
            self.k = self.j

    def step5(self):
        """step5() removes a final -e if m() > 1, and changes -ll to -l if
        m() > 1.
        """
        self.j = self.k
        if self.b[self.k] == 'e':
            a = self.m()
            if a > 1 or (a == 1 and not self.cvc(self.k-1)):
                self.k = self.k - 1
        if self.b[self.k] == 'l' and self.doublec(self.k) and self.m() > 1:
            self.k = self.k -1

    def stem(self, p, i, j):
        """In stem(p,i,j), p is a char pointer, and the string to be stemmed
        is from p[i] to p[j] inclusive. Typically i is zero and j is the
        offset to the last character of a string, (p[j+1] == '\0'). The
        stemmer adjusts the characters p[i] ... p[j] and returns the new
        end-point of the string, k. Stemming never increases word length, so
        i <= k <= j. To turn the stemmer into a module, declare 'stem' as
        extern, and delete the remainder of this file.
        """
        # copy the parameters into statics
        self.b = p
        self.k = j
        self.k0 = i
        if self.k <= self.k0 + 1:
            return self.b # --DEPARTURE--

        # With this line, strings of length 1 or 2 don't go through the
        # stemming process, although no mention is made of this in the
        # published algorithm. Remove the line to match the published
        # algorithm.

        self.step1ab()
        self.step1c()
        self.step2()
        self.step3()
        self.step4()
        self.step5()
        return self.b[self.k0:self.k+1]

## A class to handle reading a list of stopwords and providing quick lookup
class StopwordList:
    ## Constructor
    # @param fname The filename of the stopwords to read.  File should
    # contain one word per line and already be lowercase
    def __init__(self, fname=None):
        self.words = set([])
        if fname:
            self.readFile(fname)
        

    ## Read a file that contains stopwords
    # @param fname The filename of the stopwords to read.  File should
    # contain one word per line and already be lowercase            
    def readFile(self, fname):
        words = []
        try:
            f = open(fname, 'r')
            line = f.readline()
            while line:
                word = line.strip()
                if word != "":
                    words.append(word)
                line = f.readline()
            f.close()
        except Exception, e:
            msg = "Could not read file \"%s\": %s\nWill NOT be using any stopwords"
            print msg % (fname, str(e))
        self.addWords(words)


    ## Handle queries like "word in stopwords"
    # \return True only if word is in this list
    def __contains__(self, word):
        return word in self.words


    ## Same as __contains__
    # \return True only if word is in this list
    def hasWord(self, word):
        return word in self


    ## Add a list of words to the stopwords
    # @param words A set or list of words
    def addWords(self, words):
        self.words.update(words)


    ## Get the Python list that this class operates on
    # \return a Python list of the words contained in this list
    def getList(self):
        return [word for word in self.words]


    ## Return a new StopwordList composed of the words in this one and
    # another.
    # \return A new StopwordList
    def __add__(self, other):
        sl = StopwordList()
        sl.words.update(self.words, other.words)
        return sl


    ## Determine if another stopwordlist is the same as this one
    # \return True if same, False otherwise
    def __eq__(self, other):
        return self.words == other.words


    ## Determine if another stopwordlist is the same as this one
    # \return False if same, True otherwise
    def __neq__(self, other):
        return not self == other


    ## Add a list of stopwordlists together to form one superlist
    # @param stopwordLists A list of StopwordList objects
    # \return A StopwordList that is the set that contains the union of the ones given
    @classmethod
    def join(klass, stopwordLists):
        sl = StopwordList()
        words = [swl.words for swl in stopwordLists]
        sl.words.update(*words)
        return sl

## Removes from the graph the genes that are not included in a given list of genes
# @param graph GOGeneGraph or GOGenePubmedGraph from which to remove the genes
# @param geneList The genes that should be kept if they are listed in the graph
# /return graph The graph with the unwanted genes removed
def keepGenes(graph, geneList):
    for node in graph:
        remove = list()
        for gene in graph.node[node]['data'].getPropagatedGenes():
            if gene[0] not in geneList:
                remove.append(gene)
        for gene in remove:
            graph.node[node]['data'].propGenes.remove(gene)
            
        remove = list()
        for gene in graph.node[node]['data'].getGenes():
            if gene[0] not in geneList:
                remove.append(gene)
        for gene in remove:
            graph.node[node]['data'].genes.remove(gene)
    return graph
