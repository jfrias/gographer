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

genes = set([('P05091', ''), ('Q8NE62', ''), ('P47895', ''), ('P48448', ''), ('P00326', ''), ('P43353', ''), ('P08319', ''), ('P07327', '')])
nodesByGene = set(['GO:0006066', 'GO:0008152'])
if GOGenetest.getGenesByNode('GO:0006066') != genes:
    print "Error in getGenesByNode"
if GOGenetest.getNodesByGene('P48448') != nodesByGene:
    print "Error in getNodesByGene"


GOPubmedtest = GOPubmedGraph(GOGraphtest,assoc)

pubmed = set([('2347582', ''), ('7698756', ''), ('1306115', ''), ('3466164', ''), ('7828891', ''), ('8890755', '')])
propPubmed = set([('10203017', ''), ('7959792', ''), ('9473496', ''), ('7663508', ''), ('10692104', ''), ('6548414', ''), ('3466164', ''), ('10359671', ''), ('11300766', ''), ('11254442', ''), ('9000139', ''), ('10022914', ''), ('8822211', ''), ('10373482', ''), ('1695324', ''), ('1374385', ''), ('3800996', ''), ('12882974', ''), ('3755672', ''), ('9500206', ''), ('8890755', ''), ('8353497', ''), ('9733774', ''), ('9263035', ''), ('8707889', ''), ('9295054', ''), ('7698756', ''), ('8333863', ''), ('7912130', ''), ('9207799', ''), ('11597136', ''), ('9305759', ''), ('2347582', ''), ('10101297', ''), ('8359595', ''), ('11606197', ''), ('1918003', ''), ('1306115', ''), ('9268630', ''), ('9299468', ''), ('10382971', ''), ('9463486', ''), ('11179693', ''), ('9237672', ''), ('7828891', ''), ('9207021', ''), ('2110361', '')])
if GOPubmedtest.getPubMedByNode('GO:0006066') != pubmed:
    print "Error in getPubMedByNode"
for i in pubmed:
    if GOPubmedtest.getNodesByPubMed(i[0]) != set(['GO:0006066']):
        print "Error in getNodesByPubmed"
        break
if GOPubmedtest.getPropagatedPubMedByNode('GO:0008150') != propPubmed:
    print "Error in propagatePMIDs"
