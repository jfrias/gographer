#!/usr/bin/env python

"""
This is a simple tokenizer
"""
import string
import re

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
    
                
        

    