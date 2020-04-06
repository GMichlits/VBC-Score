__author__ = 'georg.michlits'

import pickle
import numpy
import math
from scipy.stats import binom

def revcomp(dna, reverse=True, complement=True):
    DNA = dna.upper()
    a = 'ACGTRYKMSWBDHVN'
    b = 'TGCAYRMKSWVHDBN'
    dict = {a[i]:b[i] for i in range (15)}
    if reverse:
        DNA = reversed(DNA)
    if complement:
        result = [dict[i] for i in DNA]
    else: result = [i for i in DNA]
    return ''.join(result)

#Species = 'Ms'
Species = 'Hm'
#Analysis = 'gene_median'
Analysis = 'sgRNAs'

guide_order_file = open('4ktwist_VBCvalidation_v2.txt','r')
Gene_d = pickle.load(open('../../Gen_properties/' + Species + '/property_sav/Gene_d.sav', 'rb'))
guide_dict = {}
Bioscore_d = pickle.load(open('../../Gen_properties/' + Species + '/property_sav/Jun_Bioscore_3Nov_d.sav', 'rb'))


gene_bioscore_lists_gn = {}
for pos in Bioscore_d:
    gene = Gene_d[pos]
    if not gene in gene_bioscore_lists_gn:
        gene_bioscore_lists_gn[gene] = []
    gene_bioscore_lists_gn[gene].append(Bioscore_d[pos])

Bioscore_write_file = open('Hmgene_Bioscore.txt','w')
Bioscore_write_file.write('gene\tgene_average_Bioscore\tgene_max_bioscore')
gene_bioscore_gn = {}
for gene in gene_bioscore_lists_gn:
    gene_av_bioscore = numpy.mean(gene_bioscore_lists_gn[gene])
    gene_bioscore_gn[gene] = gene_av_bioscore
    Bioscore_write_file.write('\n' + gene + '\t' + str(gene_av_bioscore) + '\t' + str(max(gene_bioscore_lists_gn[gene])))