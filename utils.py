import urllib, sys, time, commands, os

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
        raise RuntimeWarning, "There were %i extra pmids downloaded: %s" % (len(extrafound), ",".join(extrafound))
        #print "There were %i extra pmids downloaded: %s" % (len(extrafound), ",".join(extrafound))
        result = False
    if len(notfound) > 0:
        raise RuntimeWarning, "There were %i pmids not downloaded: %s" % (len(notfound), ",".join(notfound))                
        #print "There were %i pmids not downloaded: %s" % (len(notfound), ",".join(notfound))                
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
                print 
                raise RuntimeError, "PMID check verification failed - some documents not downloaded"
            index += blocksize
    except Exception, e:
        raise RuntimeError, "Problem fetching all PMIDs."
