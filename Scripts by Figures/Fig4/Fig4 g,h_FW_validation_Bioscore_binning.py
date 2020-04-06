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

Species = 'Ms'
#Species = 'Hm'
#Analysis = 'gene_median'
Analysis = 'sgRNAs'

guide_order_file = open('4ktwist_VBCvalidation_v2.txt','r')
guide_Species_results = open(Species + '_results_55nt.txt','r')
Gene_d = pickle.load(open('../../Gen_properties/' + Species + '/property_sav/Gene_d.sav', 'rb'))
guide_dict = {}
Bioscore_d = pickle.load(open('../../Gen_properties/' + Species + '/property_sav/Jun_Bioscore_3Nov_d.sav', 'rb'))


gene_bioscore_lists_gn = {}
for pos in Bioscore_d:
    gene = Gene_d[pos]
    if not gene in gene_bioscore_lists_gn:
        gene_bioscore_lists_gn[gene] = []
    gene_bioscore_lists_gn[gene].append(Bioscore_d[pos])

Bioscore_write_file = open('Msgene_Bioscore.txt','w')
Bioscore_write_file.write('gene\tgene_average_Bioscore\tgene_max_bioscore')
gene_bioscore_gn = {}
for gene in gene_bioscore_lists_gn:
    gene_av_bioscore = numpy.mean(gene_bioscore_lists_gn[gene])
    gene_bioscore_gn[gene] = gene_av_bioscore
    Bioscore_write_file.write('\n' + gene + '\t' + str(gene_av_bioscore) + '\t' + str(max(gene_bioscore_lists_gn[gene])))

#### read in conservation scores per gene (calculated else where - calculations_Jan2019 GOterm analysis)
gene_conservation_file = open('Msgene_conservation.txt','r')
gene_to_genecons_gd = {}
for i,line in enumerate(gene_conservation_file):
    if i > 0:
        col = line.rstrip('\n').split('\t')
        genename = col[0]
        gene_cons = float(col[1])
        gene_to_genecons_gd[genename] = gene_cons

top20VBC_d = {}
top20D16_d = {}
VBClib_d = {}
Wanglib_d = {}
sgRNA_to_pos_d = {}
nEG_pos_d = {}
CEG_pos_d = {}
Olfr_pos_d = {}
mere_sgRNAs_to_pos_d = {}
Ms_guides_details = open(Species + '_sgRNA_4order.txt','r')
for i,line in enumerate(Ms_guides_details):
    if i > 0 and not line == '\n':
        column = line.rstrip('\n').split('\t')
        if Species == 'Ms':
            genename = column[0]
            pos = column[1]
            nt30seq = column[2]
            sgRNAorder = column[3]
            purpose = column[4]
            group = column[5]
        if Species == 'Hm':
            entrezID = column[0]
            genename = column[1]
            pos = column[2]
            nt30seq = column[3]
            sgRNAorder = column[4]
            purpose = column[5]
            group = column[6]
        #print(group)
        ##### mere is not present in this lits get from excel sheet
        if 'ctrl_neg_nEG' in purpose or 'neg_ctrl' in purpose:
            nEG_pos_d[pos] = 0
        if 'ctrl_neg_olfr' in purpose:
            Olfr_pos_d[pos] = 0
        if 'ctr_pos_CEG' in purpose or 'pos_ctrl' in purpose:
            CEG_pos_d[pos] = 0
        if 'mere' in purpose:
            print('mere found')
            mere_sgRNAs_to_pos_d[sgRNAorder] = pos
        if not 'mere' in purpose:
            sgRNA_to_pos_d[sgRNAorder] = pos
            #print(pos + '\t' + sgRNAorder)
            if 'top20VBC' in group and 'top20D16' in group:
                #present in both top 20 of D16 and VBC, deactivate if you want to look only at differential sgRNAs
                top20VBC_d[pos] = 0
                top20D16_d[pos] = 0
            elif 'top20VBC' in group:
                top20VBC_d[pos] = 0
            elif 'top20D16' in group:
                top20D16_d[pos] = 0
            elif 'VBClib' in group:
                VBClib_d[pos] = 0
            elif 'Wang' in group:
                Wanglib_d[pos] = 0
            else:
                print(line)
                print('else should never be true - else CHECK CODE OR INPUT')

guide_Species_results = open(Species + '_results_55nt.txt','r')
total_reads = {}
total_reads['plasmid'] = 0
total_reads['repl1'] = 0
total_reads['repl2'] = 0
total_reads['repl3'] = 0
for i,line in enumerate(guide_Species_results):
    if i > 0 and not line == '\n':
        column = line.rstrip('\n').split('\t')
        if Species == 'Ms':
            total_reads['plasmid'] += int(column[9])
            total_reads['repl1'] += int(column[10])
            total_reads['repl2'] += int(column[11])
            total_reads['repl3'] += int(column[12])
        if Species == 'Hm':
            total_reads['plasmid'] += int(column[3])
            total_reads['repl1'] += int(column[4])
            total_reads['repl2'] += int(column[5])
            total_reads['repl3'] += int(column[6])

print(total_reads)

Structure_d = {}
#Structure_d[gene][guide1[pos,reads,LFCabs], guide2, guide3,....][guide1,guide2,guide3,guide4...]
guide_Species_results.seek(0)
for i,line in enumerate(guide_Species_results):
    if i > 0 and not line == '\n':
        column = line.rstrip('\n').split('\t')
        if Species == 'Ms':
            ID_new = column[6]
            Assigned_IDs = int(column[8])
            plas_reads = int(column[9])
            repl1_reads = int(column[10])
            repl2_reads = int(column[11])
            repl3_reads = int(column[12])
        if Species == 'Hm':
            ID_new = column[0]
            plas_reads = int(column[3])
            repl1_reads = int(column[4])
            repl2_reads = int(column[5])
            repl3_reads = int(column[6])
        plas_RP10M = (plas_reads/total_reads['plasmid']*10000000)+0.5
        repl1_RP10M = (repl1_reads/total_reads['repl1']*10000000)+0.5
        repl2_RP10M = (repl1_reads/total_reads['repl2']*10000000)+0.5
        repl3_RP10M = (repl1_reads/total_reads['repl3']*10000000)+0.5
        RP10M = float(numpy.median([repl1_RP10M,repl2_RP10M,repl3_RP10M]))
        LFC = math.log((RP10M/plas_RP10M),2)
        pvalue = binom.cdf(numpy.median([repl1_reads,repl2_reads,repl3_reads])+10,numpy.median([total_reads['repl1'],
                            total_reads['repl2'],total_reads['repl3']]),(plas_reads+10)/total_reads['plasmid'])
        if pvalue == 0.0:
            pvalue = float(1e-300)
        #print(line + '\t' + str(LFC) + '\t' + str(pvalue))
        #print(str(LFC) + '\t' + str(pvalue))
        #print(ID_new + '\t' + str(LFC))

        genename = ID_new.split('_')[0]
        sgRNA = ID_new.split('_')[1]
        if sgRNA in sgRNA_to_pos_d and Assigned_IDs == 1 and plas_reads >= 100:
            pos = sgRNA_to_pos_d[sgRNA]
            gene = Gene_d[pos]
            if Species == 'Ms':
                if gene not in Structure_d:
                    Structure_d[gene] = [[],[]]
            if Species == 'Hm':
                if gene not in Structure_d:
                    Structure_d[gene] = [[],[],[],[]]
            if pos in top20VBC_d:
                a=0
                Structure_d[gene][0].append([pos,numpy.median([repl1_reads,repl2_reads,repl3_reads]),LFC,pvalue])
                #print(genename + '\t' + str(len(sgRNA)) +'\t'+ 'top20VBC')
            if pos in top20D16_d:
                a=0
                Structure_d[gene][1].append([pos,numpy.median([repl1_reads,repl2_reads,repl3_reads]),LFC,pvalue])
                #print(genename + '\t' + str(len(sgRNA)) +'\t'+ 'top20D16')
            if pos in VBClib_d:
                a=0
                Structure_d[gene][2].append([pos,numpy.median([repl1_reads,repl2_reads,repl3_reads]),LFC,pvalue])
            if pos in Wanglib_d:
                a=0
                Structure_d[gene][3].append([pos,numpy.median([repl1_reads,repl2_reads,repl3_reads]),LFC,pvalue])

for i,gene in enumerate(Structure_d):
    if len(Structure_d[gene][0]) >= 2 and len(Structure_d[gene][1]) >= 2:
        pos = Structure_d[gene][0][0][0]
        if pos in nEG_pos_d:
            group = 'nEG'
        elif pos in CEG_pos_d:
            group = 'CEG'
        elif pos in Olfr_pos_d:
            group = 'Olfr'
        else:
            group = 'test'
        #print(str(i) +'\t' + gene + '\t' + group)
        LFC_list = []
        for element in Structure_d[gene][0]:
            pos = element[0]
            LFC = element[2]
            LFC_list.append(LFC)
        medianVBC20_LFC = numpy.median(LFC_list)
        LFCd16_list = []
        for element in Structure_d[gene][1]:
            pos = element[0]
            LFC = element[2]
            LFCd16_list.append(LFC)
        mediand16_LFC = numpy.median(LFCd16_list)
        LFCVBClib_list = []
        if Species == 'Hm':
            for element in Structure_d[gene][2]:
                pos = element[0]
                LFC = element[2]
                LFCVBClib_list.append(LFC)
            mediandVBClib_LFC = numpy.median(LFCVBClib_list)
            LFCWanglib_list = []
            for element in Structure_d[gene][1]:
                pos = element[0]
                LFC = element[2]
                LFCWanglib_list.append(LFC)
            mediandWanglib_LFC = numpy.median(LFCWanglib_list)
        if Species == 'Ms' and Analysis == 'gene_median':
            print(str(i) + '\t' + gene + '\t' + group + '\t' + str(medianVBC20_LFC) + '\t' + str(mediand16_LFC))
        if Species == 'Ms' and Analysis == 'sgRNAs':
            print(str(i) + '\t' + gene + '\t' + group + '_best' + '\t' + str(min(LFC_list)) + '\t' + str(min(LFCd16_list)) +
                  '\t' + str(gene_to_genecons_gd[gene]) + '\t' + str(numpy.median([medianVBC20_LFC,mediand16_LFC]))
                  + '\t' + str(gene_bioscore_gn[gene]))
            print(str(i) + '\t' + gene + '\t' + group + '_worst' + '\t' + str(max(LFC_list)) + '\t' + str(max(LFCd16_list)) +
                  '\t' + str(gene_to_genecons_gd[gene]) + '\t' + str(numpy.median([medianVBC20_LFC,mediand16_LFC]))
                  + '\t' + str(gene_bioscore_gn[gene]))
        if Species == 'Hm':
            print(str(i) + '\t' + gene + '\t' + group + '\t' + str(medianVBC20_LFC) + '\t' + str(mediand16_LFC)
                  + '\t' + str(mediandVBClib_LFC) + '\t' + str(mediandWanglib_LFC))
pickle.dump(Structure_d, open(Species + '_Valid_datadict.sav', 'wb'))

####
#DUMP STRUCTURE DICT
#Structure_d[gene][LIST_VBCtop20, LIST_D16top20, LIST_VBClib, LIST_Wanglib]
###  LIST VBC top:
###  [guide1[pos,reads,LFCabs,pvalue_depl], guide2, guide3,....][guide1,guide2,guide3,guide4...]



'''
############
Copy and past printed results for excel sort plot in prism
#####
'''

##########
#
#   Determine differential LFC of VBCtop2 sgRNAs versus Doench2016top2 sgRNAs)
#
##########
'''
if Species == 'Ms':
    for i,gene in enumerate(Structure_d):
        if len(Structure_d[gene][0]) >=2 and len(Structure_d[gene][1]) >= 2:
            pos = Structure_d[gene][0][0][0]
            if pos in nEG_pos_d:
                group = 'nEG'
            elif pos in CEG_pos_d:
                group = 'CEG'
            elif pos in Olfr_pos_d:
                group = 'Olfr'
            else:
                group = 'test'
            #print(str(i) +'\t' + gene + '\t' + group)
            LFC_list = []
            for element in Structure_d[gene][0]:
                pos = element[0]
                LFC = element[2]
                LFC_list.append(LFC)
            medianVBC20_LFC = numpy.median(LFC_list)
            LFCd16_list = []
            for element in Structure_d[gene][1]:
                pos = element[0]
                LFC = element[2]
                LFCd16_list.append(LFC)
            mediand16_LFC = numpy.median(LFCd16_list)
            median_all_LFC = numpy.median([mediand16_LFC, medianVBC20_LFC])
            print(str(i) + '\t' + gene + '\t' + group + '\t' + str(medianVBC20_LFC) + '\t' + str(mediand16_LFC) + '\t')

if Species == 'Hm':
    for i,gene in enumerate(Structure_d):
        if len(Structure_d[gene][0]) >=2 and len(Structure_d[gene][1]) >= 2:
            pos = Structure_d[gene][0][0][0]
            if pos in nEG_pos_d:
                group = 'nEG'
            elif pos in CEG_pos_d:
                group = 'CEG'
            elif pos in Olfr_pos_d:
                group = 'Olfr'
            else:
                group = 'test'
        print(gene + '\t' + gene)
'''

