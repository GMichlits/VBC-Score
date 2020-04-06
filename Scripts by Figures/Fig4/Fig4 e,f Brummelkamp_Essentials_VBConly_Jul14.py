__author__ = 'georg.michlits'

import matplotlib.pyplot as plt
import numpy as np
import pickle

#KBM7_Wang2015 = '../../Gen_data/v3_DICTIONARIES/Wang2015_Exp1_KBM7_abs_LFC_d.sav'
#KBM7_Wang2015 = '../../Gen_data/v3_DICTIONARIES/Wang2014_Exp17_KBM7_abs_LFC_d.sav'
KBM7_VBC = '../../Gen_data/Zuber_v3/KBM7_abs_LFC_d.sav'
KBM7_Wang2015 = '../../Gen_data/empty_data_d.sav'
#KBM7_VBC = '../../Gen_data/Zuber_v3/KBM7_abs_LFC_d.sav'
#KBM7_Wang2015 = '../../Gen_data/v3_DICTIONARIES/Wang2014_Exp17_KBM7_abs_LFC_d.sav'

print('load copynumbers')
CN_d = pickle.load(open('../../Gen_properties/Hm/property_sav/CN_d.sav', 'rb'))
print('load gene book')
Gene_d = pickle.load(open('../../Gen_properties/Hm/property_sav/Gene_d.sav', 'rb'))
print('load CEG_hart')
CEG_d = pickle.load(open('../../Gen_properties/Hm/property_sav/CEG_pos_d.sav', 'rb'))
nonEG_d = pickle.load(open('../../Gen_properties/Hm/property_sav/nonEG_d.sav', 'rb'))

print('load VBC_scores')
#VBC_score_d = pickle.load(open('../../Gen_properties/Hm/property_sav/Apr8_VBC_score_all_rel_v2_d.sav', 'rb'))
# run different scorerers to compare but also call them VBC_score_d to fool the script
'Apr8_VBC_score_all_rel_v2_d.sav'
'Jun_VBC_score_Nov3noD16_onVBC_d.sav'
#VBC_score_d = pickle.load(open('../../Gen_properties/Hm/property_sav/', 'rb'))
#VBC_score_d = pickle.load(open('../../Gen_properties/Hm/property_sav/D16_woAA_d', 'rb'))
VBC_score_d = pickle.load(open('../../Gen_properties/Hm/property_sav/Jul14_VBC_d.sav', 'rb'))
gene_cons_gd = pickle.load(open('../../Gen_properties/Hm/property_sav/AA7Genes_gd.sav', 'rb'))
print('load gene conservation groups')
gene_cons_groups = pickle.load(open('median_gene_depl_.sav', 'rb'))

reduced_to_n_sgRNAs = 3

# x = []
# n = 0
# i = 0
# k = 0
# for pos in VBC_score_d:
#     n += 1
#     x.append(VBC_score_d[pos])
#     if VBC_score_d[pos] >= 0.46:
#         i += 1
#     if VBC_score_d[pos] >= 0.50:
#         k += 1
# print(n)
# print(i/n)
# print(k/n)
# plt.hist(x,100)
# plt.show()

def return_gene_level_dg_and_CEG(screen1,screen2,CN_d,CEG_d):

    screen1_d = pickle.load(open(screen1, 'rb'))
    CN_excluded = 0
    guides_count = 0
    gene_count = 0
    genes_included = 0
    gene_level_d_sc1 = {}
    gene_score_sc1_d = {}
    gene_to_pos = {}
    for pos in screen1_d:
        if pos in VBC_score_d:
            guides_count += 1
            LFC_abs = screen1_d[pos]
            gene = Gene_d[pos]
            gene_to_pos[gene] = pos
            CN = CN_d[pos]
            if not gene in gene_level_d_sc1:
                gene_level_d_sc1[gene] = [[],[],[]]
                gene_count += 1
            if CN == 1:
                gene_level_d_sc1[gene][0].append(LFC_abs)
                gene_level_d_sc1[gene][1].append(pos)
                gene_level_d_sc1[gene][2].append(VBC_score_d[pos])
            else:
                CN_excluded += 1
    for gene in gene_level_d_sc1:
        if len(gene_level_d_sc1[gene][0]) >= 2:
            gene_score_sc1_d[gene] = np.mean(gene_level_d_sc1[gene][0])
            genes_included += 1
    print('screen1: ' + screen1)
    print('guides: ' + str(guides_count))
    print('genes: ' + str(gene_count))
    print('guides_excluded for copy number: ' + str(CN_excluded))
    print('genes with at least 2 guides: ' + str(genes_included))

    screen2_d = pickle.load(open(screen2, 'rb'))
    CN_excluded = 0
    guides_count = 0
    gene_count = 0
    genes_included = 0
    gene_level_d_sc2 = {}
    gene_score_sc2_d = {}
    for pos in screen2_d:
        if pos in VBC_score_d:
            guides_count += 1
            LFC_abs = screen2_d[pos]
            gene = Gene_d[pos]
            CN = CN_d[pos]
            if not gene in gene_level_d_sc2:
                gene_level_d_sc2[gene] = [[],[],[]]
                gene_count += 1
            if CN == 1:
                gene_level_d_sc2[gene][0].append(LFC_abs)
                gene_level_d_sc2[gene][1].append(pos)
                gene_level_d_sc2[gene][2].append(VBC_score_d[pos])
            else:
                CN_excluded += 1
    for gene in gene_level_d_sc2:
        if len(gene_level_d_sc2[gene][0]) >= 2:
            gene_score_sc2_d[gene] = np.mean(gene_level_d_sc2[gene][0])
            genes_included += 1
    print('screen2: ' + screen2)
    print('guides: ' + str(guides_count))
    print('genes: ' + str(gene_count))
    print('guides_excluded for copy number: ' + str(CN_excluded))
    print('genes with at least 2 guides: ' + str(genes_included))

    x_CEG = []
    y_CEG = []
    x_nonEG = []
    y_nonEG = []
    x = []
    y = []
    text_out_file = open('KBM7VBCvsKBM7Wang2015.txt','w')
    #text_out_file.write('gene' + '\t' + 'LFC_' +screen1 + '\t' + 'LFC_' +screen2 + '\t' + 'group')
    for gene in gene_score_sc1_d:
        #if gene in gene_score_sc2_d:
        #text_out_file.write('\n' + gene + '\t' + str(gene_score_sc1_d[gene]) + '\t' + str(gene_score_sc2_d[gene]))
        if gene_to_pos[gene] in nonEG_d:
            #x_nonEG.append(gene_score_sc2_d[gene])
            y_nonEG.append(gene_score_sc1_d[gene])
            text_out_file.write('\t' + 'non_essential')
        elif gene_to_pos[gene] in CEG_d:
            #x_CEG.append(gene_score_sc2_d[gene])
            y_CEG.append(gene_score_sc1_d[gene])
            text_out_file.write('\t' + 'core_essential')
        else:
            #x.append(gene_score_sc2_d[gene])
            y.append(gene_score_sc1_d[gene])
            text_out_file.write('\t' + 'no_group')
    #CEG2 = np.median(x_CEG)
    CEG1 = np.median(y_CEG)
    return gene_level_d_sc1, gene_level_d_sc2, CEG1

gene_level_d_sc1, gene_level_d_sc2, CEG1 = return_gene_level_dg_and_CEG(KBM7_VBC,KBM7_Wang2015,CN_d,CEG_d)

print('CEG1' + str(CEG1))

def scale_depletion(gene_level_d,CEG,CEG_reference):
    gene_level_d_new = {}
    for gene in gene_level_d:
        gene_level_d_new[gene] = gene_level_d[gene]
        for i, element in enumerate(gene_level_d[gene][0]):
            LFC_old = element
            LFC_new = LFC_old/CEG*CEG_reference
            gene_level_d_new[gene][0][i] = LFC_new
    return gene_level_d_new

#gene_level_d_sc2 = scale_depletion(gene_level_d_sc2,CEG2,CEG1)

GT_screen = open('Brummelkamp_essentials.txt','r')

# out = open('genes_mounted.txt','w')
# for gene in gene_cons_gd:
#     out.write('\n'+gene)

Brummelkamp = {}
for i, line in enumerate(GT_screen):
    if i > 0:
        col = line.rstrip('\n').split('\t')
        gene_id = col[0]
        pval = col[9]
        Brummelkamp[gene_id] = float(pval)

merged_gene_dict = {}
l = 0
k = 0
for gene in gene_level_d_sc1:
    l += 1
    #if gene in gene_level_d_sc2:
    k += 1
    merged_gene_dict[gene] = gene_level_d_sc1[gene]
    #for i,element in enumerate(gene_level_d_sc2[gene][0]):
    #    merged_gene_dict[gene][0].append(gene_level_d_sc2[gene][0][i])
    #    merged_gene_dict[gene][1].append(gene_level_d_sc2[gene][1][i])
    #    merged_gene_dict[gene][2].append(gene_level_d_sc2[gene][2][i])
print('screen1guides: ' + str(l))
#print('screen2guides: ' + str(k))

i = 0
genes_pop = []
sgRNA_per_gene_list = []
for gene in merged_gene_dict:
    sgRNA_per_gene_list.append(len(merged_gene_dict[gene][0]))
    if len(merged_gene_dict[gene][0]) < 6:
        i += 1
        genes_pop.append(gene)
for gene in genes_pop:
    merged_gene_dict.pop(gene)
print('popped genes: ' + str(i))
print('average sgRNAs per gene' + str(np.mean(sgRNA_per_gene_list)))
############## info on structure of dictionary
# gene_level_d_sc1[gene][0].append(LFC_abs)
# gene_level_d_sc1[gene][1].append(pos)
# gene_level_d_sc1[gene][2].append(VBC_score_d[pos])

top4_VBCs_gene_dict = {}
rand4_gene_dict = {}
for gene in merged_gene_dict:
    VBCs_gene = []
    for VBC_score in merged_gene_dict[gene][2]:
        VBCs_gene.append(VBC_score)
    i = 0
    top4_VBCs = []
    for score in sorted(VBCs_gene, reverse=True):
        if i < reduced_to_n_sgRNAs:
            top4_VBCs.append(score)
        i += 1
        if i == reduced_to_n_sgRNAs:
            break
    top4_VBCs_gene_dict[gene] = [[],[],[]]
    rand4_gene_dict[gene] = [[],[],[]]
    do_not_repeat_rand_index = []
    for top4_score in top4_VBCs:
        test_new_random = 0
        while test_new_random == 0:
            random_guide_index = int(np.random.random()*len(merged_gene_dict[gene][2]))
            if not random_guide_index in do_not_repeat_rand_index:
                do_not_repeat_rand_index.append(random_guide_index)
                test_new_random += 1
        top_guide_index = merged_gene_dict[gene][2].index(top4_score)
        top4_VBCs_gene_dict[gene][0].append(merged_gene_dict[gene][0][top_guide_index])
        top4_VBCs_gene_dict[gene][1].append(merged_gene_dict[gene][1][top_guide_index])
        top4_VBCs_gene_dict[gene][2].append(merged_gene_dict[gene][2][top_guide_index])
        rand4_gene_dict[gene][0].append(merged_gene_dict[gene][0][random_guide_index])
        rand4_gene_dict[gene][1].append(merged_gene_dict[gene][1][random_guide_index])
        rand4_gene_dict[gene][2].append(merged_gene_dict[gene][2][random_guide_index])

mean_LFC_dict = {}
mean_VBC_dict = {}
for gene in merged_gene_dict:
    mean_LFC = np.mean(merged_gene_dict[gene][0])
    mean_VBC = np.mean(merged_gene_dict[gene][2])
    mean_LFC_dict[gene] = mean_LFC
    mean_VBC_dict[gene] = mean_VBC

top4_mean_LFC_dict = {}
top4_mean_VBC_dict = {}
rand4_mean_LFC_dict = {}
for gene in top4_VBCs_gene_dict:
    top4_mean_LFC = np.mean(top4_VBCs_gene_dict[gene][0])
    rand4_mean_LFC = np.mean(rand4_gene_dict[gene][0])
    top4_mean_VBC = np.mean(top4_VBCs_gene_dict[gene][2])
    rand4_mean_LFC_dict[gene] = rand4_mean_LFC
    top4_mean_LFC_dict[gene] = top4_mean_LFC
    top4_mean_VBC_dict[gene] = top4_mean_VBC

#this paragraph was only used to optimize cutoff criteria
'''
pvle_cutoff_BK = [3,5,7,9,11,13,15,17,20,22,25,30,40]
for pvle_cutoff in pvle_cutoff_BK:
    i = 0
    for gene in Brummelkamp:
        if Brummelkamp[gene] >= pvle_cutoff:
            i+=1
    # print('pvalue_cutoff: ' + str(pvle_cutoff_BK) + '\t' 'essentials: ' + str(i))

cutoff_essential_in_CR = [-1,CEG1/3,-1.5,-1.7,-2,-2.5,-3]
for cutoff in cutoff_essential_in_CR:
    i = 0
    for gene in mean_LFC_dict:
        if mean_LFC_dict[gene] <= cutoff:
            i += 1
    # print('cutoff: ' + str(cutoff) + '\t' + 'essential genes: ' + str(i))
'''
def missed_overlap_allsgRNA(GT_screens,CRISPR_screens,GT_cutoff,CR_cutoff):
    BK_ess = []
    for gene in GT_screens:
        if gene in top4_mean_LFC_dict:
            if GT_screens[gene] >= GT_cutoff:
                BK_ess.append(gene)

    CR_ess = []
    for gene in CRISPR_screens:
        if gene in top4_mean_LFC_dict:
            if CRISPR_screens[gene] <= CR_cutoff:
                CR_ess.append(gene)
    overlap = []
    overlap2 = []
    only_CR = []
    only_GT = []
    for gene in BK_ess:
        if gene in top4_mean_LFC_dict:
            if gene in CR_ess:
                overlap.append(gene)
            else:
                only_GT.append(gene)
    for gene in CR_ess:
        if gene in top4_mean_LFC_dict:
            if gene in BK_ess:
                overlap2.append(gene)
            else:
                only_CR.append(gene)
    print('GeneTrap essentials: ' + str(len(BK_ess)))
    print('CRISPR essentials: ' + str(len(CR_ess)))
    print('Only_GT: ' + str(len(only_GT)))
    print('Only_CR: ' + str(len(only_CR)))
    print('overlap: ' + str(len(overlap)))
    print('overlap2_confirmation: ' + str(len(overlap2)))
    return BK_ess, CR_ess, only_GT, only_CR, overlap
def missed_overlap_top4VBC(GT_screens,CRISPR_screens,GT_cutoff,CR_cutoff):
    BK_ess = []
    for gene in GT_screens:
        if gene in top4_mean_LFC_dict:
            if GT_screens[gene] >= GT_cutoff:
                BK_ess.append(gene)

    CR_ess = []
    for gene in CRISPR_screens:
        if gene in top4_mean_LFC_dict:
            if top4_mean_LFC_dict[gene] <= CR_cutoff:
                CR_ess.append(gene)
    overlap = []
    overlap2 = []
    only_CR = []
    only_GT = []
    for gene in BK_ess:
        if gene in top4_mean_LFC_dict:
            if gene in CR_ess:
                overlap.append(gene)
            else:
                only_GT.append(gene)
    for gene in CR_ess:
        if gene in top4_mean_LFC_dict:
            if gene in BK_ess:
                overlap2.append(gene)
            else:
                only_CR.append(gene)
    print('GeneTrap essentials: ' + str(len(BK_ess)))
    print('CRISPR essentials: ' + str(len(CR_ess)))
    print('Only_GT: ' + str(len(only_GT)))
    print('Only_CR: ' + str(len(only_CR)))
    print('overlap: ' + str(len(overlap)))
    print('overlap2_confirmation: ' + str(len(overlap2)))
    return BK_ess, CR_ess, only_GT, only_CR, overlap
def missed_overlap_rand4(GT_screens,CRISPR_screens,GT_cutoff,CR_cutoff):
    BK_ess = []
    for gene in GT_screens:
        if gene in top4_mean_LFC_dict:
            if GT_screens[gene] >= GT_cutoff:
                BK_ess.append(gene)

    CR_ess = []
    for gene in CRISPR_screens:
        if gene in top4_mean_LFC_dict:
            if rand4_mean_LFC_dict[gene] <= CR_cutoff:
                CR_ess.append(gene)
    overlap = []
    overlap2 = []
    only_CR = []
    only_GT = []
    for gene in BK_ess:
        if gene in top4_mean_LFC_dict:
            if gene in CR_ess:
                overlap.append(gene)
            else:
                only_GT.append(gene)
    for gene in CR_ess:
        if gene in top4_mean_LFC_dict:
            if gene in BK_ess:
                overlap2.append(gene)
            else:
                only_CR.append(gene)
    print('GeneTrap essentials: ' + str(len(BK_ess)))
    print('CRISPR essentials: ' + str(len(CR_ess)))
    print('Only_GT: ' + str(len(only_GT)))
    print('Only_CR: ' + str(len(only_CR)))
    print('overlap: ' + str(len(overlap)))
    print('overlap2_confirmation: ' + str(len(overlap2)))
    return BK_ess, CR_ess, only_GT, only_CR, overlap


# for GT_c in pvle_cutoff_BK:
#     for LFC_cut in cutoff_essential_in_CR:
GT_c = 9
LFC_cut = CEG1/3

print('TOP_4_VBC_sgRNA')
BK_ess, CR_ess, only_GT, only_CR, overlap = missed_overlap_top4VBC(Brummelkamp, mean_LFC_dict,GT_c,LFC_cut)
print('rand4_VBC_sgRNA')
BK_ess, CR_ess, only_GT, only_CR, overlap = missed_overlap_rand4(Brummelkamp, mean_LFC_dict,GT_c,LFC_cut)
print('ALL_16_sgRNA')
BK_ess, CR_ess, only_GT, only_CR, overlap = missed_overlap_allsgRNA(Brummelkamp, mean_LFC_dict,GT_c,LFC_cut)


print('BK_ess')
print(BK_ess)
print('CR_ess')
print(CR_ess)

'''
# print('only_GT')
# print(only_GT)
# print('only_CR')
# print(only_CR)
# print('overlap')
# print(overlap)
groups = [only_GT, only_CR, overlap]
for group in groups:
    VBC_scores = []
    for gene in group:
        if gene in merged_gene_dict:
            VBC_scores.append(np.mean(merged_gene_dict[gene][2]))
    average_VBC_scores = np.mean(VBC_scores)
    # print(average_VBC_scores)

i = 0
k = 0
l = 0
lo = 0
me = 0
hi = 0

#commented out for conservation variance
'''
onlyGT_cons = {}
for gene in only_GT:
    if gene in gene_cons_gd:
        cons = gene_cons_gd[gene]
        onlyGT_cons[gene] = cons

only_GT_hicons = []
only_GT_mecons = []
only_GT_locons = []
total_len = len(onlyGT_cons)
hi_med = total_len/3
med_low = total_len/3*2
for i, item in enumerate(sorted(onlyGT_cons.items(), key=lambda t: t[1], reverse = True)):
    gene = item[0]
    cons = item[1]
    if i < hi_med:
        only_GT_hicons.append(gene)
    elif i < med_low:
        only_GT_mecons.append(gene)
    else:
        only_GT_locons.append(gene)
'''
onlyGT_deltaVBC = {}
for gene in only_GT:
    if gene in top4_mean_VBC_dict:
        delta_VBC = top4_mean_VBC_dict[gene]/mean_VBC_dict[gene]
        onlyGT_deltaVBC[gene] = delta_VBC

only_GT_hicons = []
only_GT_mecons = []
only_GT_locons = []
total_len = len(onlyGT_deltaVBC)
hi_med = total_len/3
med_low = total_len/3*2
VBC_score_list = []
for i, item in enumerate(sorted(onlyGT_deltaVBC.items(), key=lambda t: t[1], reverse = True)):
    gene = item[0]
    cons = item[1]
    VBC_score_list.append(cons)
    if i < hi_med:
        only_GT_hicons.append(gene)
    elif i < med_low:
        only_GT_mecons.append(gene)
    else:
        only_GT_locons.append(gene)

i = 0
for gene in only_GT:
    #if gene in mean_LFC_dict and gene in top4_mean_LFC_dict:
    i += 1
    if top4_mean_LFC_dict[gene] < LFC_cut:
        k += 1
        if gene in only_GT_hicons:
            hi += 1
        elif gene in only_GT_mecons:
            me += 1
        elif gene in only_GT_locons:
            lo += 1
    else:
        l += 1
print(i)
print(k)
print(l)
print('high, medium ,low')
print(hi)
print(me)
print(lo)

print('GT_pvalue ' + '\t' + str(GT_c) + '\t' + 'LFC_cutoff ' +'\t' + str(LFC_cut) +'\t' + 'missed genes ' +'\t' + str(i) + '\t' + 'rescure rate ' +'\t' + str(k/i) + '\t' + 'low vs hi ' + '\t' + str(lo/hi))
'''