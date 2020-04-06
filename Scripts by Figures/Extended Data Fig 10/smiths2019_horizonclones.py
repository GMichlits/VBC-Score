__author__ = 'georg.michlits'

import pickle
import numpy as np

horfile = open('horizonClones3.txt', 'r')
b1file = open('1b.txt', 'r')

Species = 'Hm'
lookup_d = pickle.load(open('../../Gen_properties/' + Species + '/property_sav/Gene_to_20nt_to_pos.sav', 'rb'))
Bioscore_d = pickle.load(open('../../Gen_properties/' + Species + '/property_sav/Jul14_Bioscore_d.sav', 'rb'))
distance_d = pickle.load(open('../../Gen_properties/' + Species + '/property_sav/distance_d.sav', 'rb'))
Gene_d = pickle.load(open('../../Gen_properties/' + Species + '/property_sav/Gene_d.sav', 'rb'))

HAP1_screen = (pickle.load(open('../../Gen_data/v3_DICTIONARIES/TKOv3_HAP1d18_abs_d.sav', 'rb')))

gene_LFC_guides = {}
for pos in HAP1_screen:
    gene = Gene_d[pos]
    if not gene in gene_LFC_guides:
        gene_LFC_guides[gene] = []
    gene_LFC_guides[gene].append(HAP1_screen[pos])

gene_LFC_d = {}
for gene in gene_LFC_guides:
    gene_LFC_d[gene] = np.mean(gene_LFC_guides[gene])

#get_RNA[clone] = RNA
get_RNA = {}
for i,line in enumerate(horfile):
    if i > 0:
        col = line.rstrip('\n').split('\t')
        if col[3] == 'RNA':
            clone = col[0]
            get_RNA[clone] = col[4]

horfile.seek(0)
k = 0
for i,line in enumerate(horfile):
    if i > 0:
        col = line.rstrip('\n').split('\t')
        if len(col) > 5:
            gene = col[1]
            clone = col[0]
            sgRNA = col[5]
            if not sgRNA == '' and not sgRNA == 'n.f.':
                if gene in lookup_d:
                    if sgRNA in lookup_d[gene]:
                        pos = lookup_d[gene][sgRNA]
                        Bioscore = Bioscore_d[pos]
                        distance = distance_d[pos]
                        if gene in gene_LFC_d:
                            LFC_gene = gene_LFC_d[gene]
                        else:
                            LFC_gene = 'na'
                        print(line.rstrip('\n') + '\t' + str(Bioscore) + '\t' + str(distance) + '\t' + str(LFC_gene) + '\t' + str(get_RNA[clone]))
                        k += 1
                    #else:
                        #print(line.rstrip('\n') + '\t' + 'sgRNA not found' + '\t' + 'sgRNA not found' + '\t' + str(LFC_gene))
                #else:
                    #print(line.rstrip('\n') + '\t' + 'gene not found' + '\t' + 'gene not found' + '\t' + 'gene not found')
print(k)

