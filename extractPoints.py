f = open('pointsIBrandProtein.txt')
points = dict()

##format of files
#Column 1: # of proteins
#Column 2: # of directly associated nodes
#Column 3: # of nodes in Steiner tree
#Column 4: Length of Steiner tree

for line in f:
    data = line.split('\t')
    if int(data[0]) not in points:
        points[int(data[0])] = list()
    points[int(data[0])].append(float(data[3]))

f = open('geneRandProcessedDict030612.txt')
pointsPPI = dict()

##format of files
#Column 1: # of genes
#Column 2: Info loss
#Column 3: Probability
#Column 4: # possible PPI
#Column 5: # found PPI

for line in f:
    data = line.strip().split('\t')
    if int(data[0]) not in pointsPPI:
        pointsPPI[int(data[0])] = list()
    pointsPPI[int(data[0])].append(float(data[4]))
