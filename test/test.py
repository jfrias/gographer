from gographer import *

location = "./testgraph.obo-xml"

#possible namespaces are: "biological_process", "molecular_function", "cellular_component"
GOGraphtest = GOGraph("biological_process", location)

descrip = 'The chemical reactions and pathways involving (R)-4-hydroxymandelate, the anion of a hydroxylated derivative of mandelate (alpha-hydroxybenzeneacetate).'
namespace = 'biological_process'
if GOGraphtest.getNodeDescription('GO:0046431') != descrip:
    print "Error in getNodeDescription"
if GOGraphtest.getNodeNamespace('GO:0046431') != namespace:
    print "Error in getNodeNamespace"

assoc = "./gene_association.test"
GOGenetest = GOGeneGraph(GOGraphtest,assoc)

genes = set([('Q8NE62', ''), ('P48448', ''), ('P43353', ''), ('P05091', ''), ('P07327', ''), ('P08319', '')])
propGenes = set([('P05091', ''), ('Q8NE62', ''), ('O43704', ''), ('P48448', ''), ('P07327', ''), ('P43353', ''), ('P08319', '')])
nodesByGene = set(['GO:0006725', 'GO:0032787'])
if GOGenetest.getGenesByNode('GO:0006066') != genes:
    print "Error in getGenesByNode"
if GOGenetest.getNodesByGene('P05177') != nodesByGene:
    print "Error in getNodesByGene"
if GOGenetest.getPropagatedGenesByNode('GO:0006066') != propGenes:
    print "Error in getPropagatedGenesByNode"


GOPubmedtest = GOPubmedGraph(GOGraphtest,assoc)

pubmed = set([('1306115', ''), ('8890755', ''), ('2347582', ''), ('3466164', ''), ('7828891', '')])
propPubmed = set([('2347582', ''), ('1306115', ''), ('9463486', ''), ('8890755', ''), ('7828891', ''), ('3466164', '')])
if GOPubmedtest.getPubMedByNode('GO:0006066') != pubmed:
    print "Error in getPubMedByNode"
for i in pubmed:
    if GOPubmedtest.getNodesByPubMed(i[0]) != set(['GO:0006066']):
        print "Error in getNodesByPubmed"
        break
if GOPubmedtest.getPropagatedPubMedByNode('GO:0006066') != propPubmed:
    print "Error in getPropagatedPubMedByNode"

corpus = Corpus.fromPubmedArticleSetFile("test_block_0.xml")
stopwords = StopwordList("stopwords.txt")
GOPubmedtest.calculateWordVectors(corpus, stopwords = stopwords)
