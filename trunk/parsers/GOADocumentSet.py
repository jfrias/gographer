from xml.sax import handler, make_parser
import re

#Remove this when making package
import sys
sys.path.append('../..')

from gographer.Document import Document

## This class parses GOA files.  When parsing is finished, the docs
# instance variable contains a dictionary of the form {'pmid': Document}
class GOADocumentSet(handler.ContentHandler):
    def __init__(self):
        handler.feature_external_ges = "false"
        self.docs = {}
        self.doc = None
        self.chars = ""
        
        
    def startElement(self, name, attr):
        if name == 'DOC':
            self.doc = Document()
            self.doc.goids= []
            self.doc.pmid = attr['PMID'].encode('ascii', 'ignore')
        self.chars = ""

            
    def endElement(self, name):
        if name == 'DOC':
            self.docs[self.doc.pmid] = self.doc
        if name == 'Abstract':
            self.doc.abstract = self.text()
        if name == "Annots":
            self.doc.goids.append(self.text())
            
    
    def characters(self, data):
        self.chars += data


    def text(self):
        return self.chars.strip().encode('ascii', 'ignore')
    

    ## Method to parse a GOADocumentSet XML file.
    # @param location The location of the xml file to parse
    # \return A GOADocumentSet object 
    @classmethod
    def parse(self, location):
        parser = make_parser()
        parser.setFeature("http://xml.org/sax/features/external-general-entities", False)
        parser.setFeature("http://xml.org/sax/features/external-parameter-entities", False)
        handler = GOADocumentSet()
        parser.setContentHandler(handler)
        try:
            f = open(location, 'r')
            parser.parse(f)
            f.close()
        except Exception, e:
            raise RuntimeError, "Could not parse GOADocumentSet XML file at %s" % location
        return handler

            
                                                                        



