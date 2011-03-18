from GOGraph import *
from GOGeneGraph import *
from GOPubmedGraph import *

location = "./go_daily-termdb.obo-xml"
location = "./testgraph.obo-xml"

#possible namespaces are: "biological_process", "molecular_function", "cellular_component"
GOGraphtest = GOGraph("biological_process", location)

descrip = 'The chemical reactions and pathways involving (R)-4-hydroxymandelate, the anion of a hydroxylated derivative of mandelate (alpha-hydroxybenzeneacetate).'
namespace = 'biological_process'
if GOGraphtest.getNodeDescription('GO:0046431') != descrip:
    print "Error in getNodeDescription"
if GOGraphtest.getNodeNamespace('GO:0046431') != namespace:
    print "Error in getNodeNamespace"

assoc = "./gene_association.goa_human"
GOGenetest = GOGeneGraph(GOGraphtest,assoc)

GOPubmedtest = GOPubmedGraph(GOGraphtest,assoc)
