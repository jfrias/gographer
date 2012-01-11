import networkx,re

m = re.search("\.",networkx.__version__)
if float(networkx.__version__[:m.start()]) < float("1.0"):
    raise RuntimeError, "You must upgrade networkx to at least version 1.0."

from utils import *
from GOGraph import *
from GOGeneGraph import *
from GOPubmedGraph import *
from GOGenePubmedGraph import *
from Corpus import *
from SteinerTree import *
