f = open('downStream.Human.txt')
geneNames = set()
for i in f:
    i = i.strip()
    i = i.split('\t')
    for gene in i:
        geneNames.add(gene)
proteinIDLookup = dict()
f = open('./Data/uniprot-organism-Homo-sapiens.tab')
for line in f:
    if line.find('Accession') == 0:
        continue
    else:
        line = line.strip()
        columns = line.split('\t')
        ID = columns[0]
        names = columns[4].split(' ')
        for name in names:
            if name not in proteinIDLookup:
                proteinIDLookup[name] = set([ID])
            else:
                proteinIDLookup[name].add(ID)
f.close()
geneIDs = set()
for i in list(geneNames):
    geneIDs = geneIDs.union(proteinIDLookup[i])

uniprotIB = GOGraph.loadPickle('GOGeneDaily.uniprot.noiea.5.19.11.IB')
uniprotIB.parseAssocFile = GOGeneGraph.parseAssocFile
uniprotIB.excludeEvidence = ['IEA']
keepGenes(uniprotIB, list())
uniprotIB.parseAssocFile(uniprotIB, 'gene_association.goa_human')
keepGenes(uniprotIB, list(geneIDs))
uniprotIB.removeGeneless()
uniprotIB.removePMIDless()
mergeGraphCheck(uniprotIB, modelIB)

uniprotIBGeneMult = GOGraph.loadPickle('GOGeneDaily.uniprot.noiea.5.19.11.IB')
uniprotIBGeneMult.parseAssocFile = GOGeneGraph.parseAssocFile
uniprotIBGeneMult.excludeEvidence = ['IEA']
keepGenes(uniprotIBGeneMult, list())
uniprotIBGeneMult.parseAssocFile(uniprotIBGeneMult, 'gene_association.goa_human')
keepGenes(uniprotIBGeneMult, list(geneIDs))
uniprotIBGeneMult.removeGeneless()
uniprotIBGeneMult.removePMIDless()
mergeGraphCheck(uniprotIBGeneMult, modelIBGeneMult)
mergeGraphMultCheck(uniprotIBGeneMult, modelIB)

uniprotIBRand = GOGraph.loadPickle('GOGeneDaily.uniprot.noiea.5.19.11.IB')
uniprotIBRand.parseAssocFile = GOGeneGraph.parseAssocFile
uniprotIBRand.excludeEvidence = ['IEA']
keepGenes(uniprotIBRand, list())
uniprotIBRand.parseAssocFile(uniprotIBRand, 'gene_association.goa_human')
keepGenes(uniprotIBRand, list(geneIDs))
uniprotIBRand.removeGeneless()
uniprotIBRand.removePMIDless()
mergeGraphCheckRand(uniprotIBRand, modelIB)

uniprotIBGeneMultRand = GOGraph.loadPickle('GOGeneDaily.uniprot.noiea.5.19.11.IB')
uniprotIBGeneMultRand.parseAssocFile = GOGeneGraph.parseAssocFile
uniprotIBGeneMultRand.excludeEvidence = ['IEA']
keepGenes(uniprotIBGeneMultRand, list())
uniprotIBGeneMultRand.parseAssocFile(uniprotIBGeneMultRand, 'gene_association.goa_human')
keepGenes(uniprotIBGeneMultRand, list(geneIDs))
uniprotIBGeneMultRand.removeGeneless()
uniprotIBGeneMultRand.removePMIDless()
mergeGraphMultCheckRand(uniprotIBGeneMultRand, modelIB)
