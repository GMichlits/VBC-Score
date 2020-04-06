__author__ = 'georg.michlits'

import numpy as np
import matplotlib.pyplot as plt
import pickle

print('load lfc screen data')
Zub_screens = []
Zub_screens.append('../../Gen_data/Zuber_v3/K0.5KBM7_abs_LFC_d.sav')
compare_screens = []
compare_screens.append([])
#compare_screens[0].append('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Wang2014_KBM7_Exp17abs_LFC_d.sav')
#compare_screens[0].append('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Exp17abs_LFC_d.sav')
#compare_screens[0].append('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Exp1abs_LFC_d.sav')
compare_screens[0].append('../../Gen_data/GenomeCrispr/Exp_data_v3/Exp1abs_LFC_d.sav')
Zub_screens.append('../../Gen_data/Zuber_v3/K0.5RKO_abs_LFC_d.sav')
compare_screens.append([])
#compare_screens[1].append('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Munoz_RKO_Exp21abs_LFC_d.sav')
#compare_screens[1].append('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Brunello_abs_LFC_d.sav')
#compare_screens[1].append('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Yusa_HT1080_abs_LFC_d.sav')
#compare_screens[1].append('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/TKOv3_HAP1d18_abs_d.sav')
compare_screens[1].append('../../Gen_data/GenomeCrispr/Exp_data_v3/Exp145abs_LFC_d.sav')
Zub_screens.append('../../Gen_data/Zuber_v3/K0.5_MIAPACA2_abs_LFC_d.sav')
compare_screens.append([])
compare_screens[2].append('../../Gen_data/GenomeCrispr/Exp_data_v3/Exp340abs_LFC_d.sav')

print('load copynumbers')
CN_d = pickle.load(open('../../Gen_properties/Hm/property_sav/CN_d.sav', 'rb'))
print('load gene book')
Gene_d = pickle.load(open('../../Gen_properties/Hm/property_sav/Gene_d.sav', 'rb'))
print('load CEG_hart')
CEG_d = pickle.load(open('../../Gen_properties/Hm/property_sav/CEG_pos_d.sav', 'rb'))
nonEG_d = pickle.load(open('../../Gen_properties/Hm/property_sav/nonEG_d.sav', 'rb'))

def savefig_screen1vs2(screen1,screen2,CN_d,CEG_d):
    screen1_d = pickle.load(open(screen1, 'rb'))
    CN_excluded = 0
    guides_count = 0
    gene_count = 0
    genes_included = 0
    gene_level_d = {}
    gene_score_sc1_d = {}
    gene_to_pos = {}
    for pos in screen1_d:
        guides_count += 1
        LFC_abs = screen1_d[pos]
        gene = Gene_d[pos]
        gene_to_pos[gene] = pos
        CN = CN_d[pos]
        if not gene in gene_level_d:
            gene_level_d[gene] = []
            gene_count += 1
        if CN == 1:
            gene_level_d[gene].append(LFC_abs)
        else:
            CN_excluded += 1
    for gene in gene_level_d:
        if len(gene_level_d[gene]) >= 3:
            gene_score_sc1_d[gene] = np.mean(gene_level_d[gene])
            genes_included += 1
    print('screen1: ' + screen1)
    print('guides: ' + str(guides_count))
    print('genes: ' + str(gene_count))
    print('guides_excluded for copy number: ' + str(CN_excluded))
    print('genes with at least 3 guides: ' + str(genes_included))

    screen2_d = pickle.load(open(screen2, 'rb'))
    CN_excluded = 0
    guides_count = 0
    gene_count = 0
    genes_included = 0
    gene_level_d = {}
    gene_score_sc2_d = {}
    for pos in screen2_d:
        guides_count += 1
        LFC_abs = screen2_d[pos]
        gene = Gene_d[pos]
        CN = CN_d[pos]
        if not gene in gene_level_d:
            gene_level_d[gene] = []
            gene_count += 1
        if CN == 1:
            gene_level_d[gene].append(LFC_abs)
        else:
            CN_excluded += 1
    for gene in gene_level_d:
        if len(gene_level_d[gene]) >= 3:
            gene_score_sc2_d[gene] = np.mean(gene_level_d[gene])
            genes_included += 1
    print('screen2: ' + screen2)
    print('guides: ' + str(guides_count))
    print('genes: ' + str(gene_count))
    print('guides_excluded for copy number: ' + str(CN_excluded))
    print('genes with at least 3 guides: ' + str(genes_included))

    x_CEG = []
    y_CEG = []
    x_nonEG = []
    y_nonEG = []
    x = []
    y = []
    text_out_file = open(screen1.split('/')[-1] + 'vs' + screen2.split('/')[-1] + '.txt','w')
    text_out_file.write('gene' + '\t' + 'LFC_' +screen1 + '\t' + 'LFC_' +screen2 + '\t' + 'group')
    for gene in gene_score_sc1_d:
        if gene in gene_score_sc2_d:
            text_out_file.write('\n' + gene + '\t' + str(gene_score_sc1_d[gene]) + '\t' + str(gene_score_sc2_d[gene]))
            if gene_to_pos[gene] in nonEG_d:
                x_nonEG.append(gene_score_sc2_d[gene])
                y_nonEG.append(gene_score_sc1_d[gene])
                text_out_file.write('\t' + 'non_essential')
            elif gene_to_pos[gene] in CEG_d:
                x_CEG.append(gene_score_sc2_d[gene])
                y_CEG.append(gene_score_sc1_d[gene])
                text_out_file.write('\t' + 'core_essential')
            else:
                x.append(gene_score_sc2_d[gene])
                y.append(gene_score_sc1_d[gene])
                text_out_file.write('\t' + 'no_group')

    CEG1 = np.median(x_CEG)
    CEG2 = np.median(y_CEG)
    print('other' + str(CEG1) + ' IMP_lib' + str(CEG2))
    plt.figure(figsize=(12,16.0))
    plt.tick_params(labelsize=40)
    plt.scatter(x,y,c='k', s=100, alpha=0.25, linewidths=None)
    #plt.plot([-8,CEG2],[4,CEG2],'r--', lw=4)
    #plt.plot([CEG1,-12],[CEG1,2.5],'r--', lw=4)
    plt.scatter(x_CEG,y_CEG,c='red', s=100, alpha=1)
    plt.scatter(x_nonEG,y_nonEG,c='b', s=100, alpha=1)
    plt.xlim(-6,4)
    plt.ylim(-10,2.5)
    plt.grid()
    plt.savefig(screen1.split('/')[-1] + '_vs_' + screen2.split('/')[-1] + '.png')
    plt.savefig(screen1.split('/')[-1] + '_vs_' + screen2.split('/')[-1] + '.pdf')

for i in range(len(Zub_screens)):
    for screen in compare_screens[i]:
        screen1 = Zub_screens[i]
        screen2 = screen
        savefig_screen1vs2(screen1,screen2,CN_d,CEG_d)
