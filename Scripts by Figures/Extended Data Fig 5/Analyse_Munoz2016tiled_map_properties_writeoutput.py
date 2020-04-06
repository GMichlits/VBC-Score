__author__ = 'georg.michlits'

import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split


print('load dataset')
#RES16_LFC_d = pickle.load(open('../../Gen_data/Novartis_tiled/DLD_abs_LFC_d.sav', 'rb'))
#RES16_LFC_d = pickle.load(open('../../Gen_data/Novartis_tiled/RKO_abs_LFC_d.sav', 'rb'))
RES16_LFC_d = pickle.load(open('../../Gen_data/Novartis_tiled/NCIH1299_abs_LFC_d.sav', 'rb'))

Species = 'Hm'
print('load Features')
distance_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/distance_d.sav', 'rb'))
VBC_score_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/D16_woAA_d.sav', 'rb'))
VBC_score_new_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/Jul14_VBC_d.sav', 'rb'))
inDelphi_infr_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/inDelphi_infr_d.sav', 'rb'))

print('load other useful dictionaries')
Gene_d = pickle.load(open('../../Gen_properties/' + Species + '/property_sav/Gene_d.sav', 'rb'))

##################### dictionaries for printed text file
nt30_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/nt30_d.sav', 'rb'))
phylo21_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/phylo21_d.sav', 'rb'))
phast21_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/phast21_d.sav', 'rb'))
Pfam_family_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/Pfam_family_d.sav', 'rb'))
Pfam_domain_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/Pfam_domain_d.sav', 'rb'))

VBC_score_rev_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/Jul14_VBC_d.sav', 'rb'))
Bioscore_rev_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/Jul14_Bioscore_d.sav', 'rb'))
tracrv2_rev_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/tracrv2_3Nov_d.sav', 'rb'))
D16_woAA_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/D16_woAA_d.sav', 'rb'))
inDelphi = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/inDelphi_infr_d.sav', 'rb'))
AA_cons_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/AA7_d.sav', 'rb'))

print('load AA_scores')
AA_score_dm2 = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/AA_residue_d-2.sav', 'rb'))
AA_score_dm1 = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/AA_residue_d-1.sav', 'rb'))
AA_score_dv2 = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/AA_residue_dv2.sav', 'rb'))
AA_score_dp1 = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/AA_residue_d+1.sav', 'rb'))
AA_score_dp2 = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/AA_residue_d+2.sav', 'rb'))
print('load scoretables')
AA_score_m2 = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/AA_score_m2.sav', 'rb'))
AA_score_m1 = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/AA_score_m1.sav', 'rb'))
AA_score_0 = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/AA_score_0.sav', 'rb'))
AA_score_p1 = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/AA_score_p1.sav', 'rb'))
AA_score_p2 = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/AA_score_p2.sav', 'rb'))
###################
#
#  NOTES
#   All dictionaries with _d contain "NAME_Property"_d = {pos : score}
#   All scores are between 0 and 1
#   pos is a guide ID string: chr12:85822107-85822137(+)
#   Dictionaries AA_cons_d has option of choosing window size "NAME_Property"_d = {pos : {1: score,
#                                                                                         3: score,
#                                                                                         5: score,.....}
####################

Data_list_dict = {}
Structure_dict = {}
outfile_tiled = open('outfile_tiled_NCIH1299.txt', 'w')
outfile_tiled.write('gene\tsgRNA_pos\tnt30seq\tabs_depletion\tdistance\tPfam_dom\tPfam_fam\tAA_cons\tphylo\tphast' +
                    '\tAAidentity\ttracrv2score\tDoench2016score\tframeshifts\tBioscore\tVBCscore')

print('generatefeatures for model')
i = 0
guide_not_found = 0
data = []
for pos in RES16_LFC_d:
    if pos in VBC_score_d:
        data_list = []
        #####################################################################################
        ##                    x Features for MODEL
        AA_cons = VBC_score_new_d[pos]
        data_list.append(AA_cons)
        data_list.append(VBC_score_d[pos])
        gene = Gene_d[pos]
        Data_list_dict[pos] = data_list
        outfile_tiled.write('\n' + gene)
        outfile_tiled.write('\t' + pos)
        outfile_tiled.write('\t' + nt30_d[pos])
        outfile_tiled.write('\t' + str(RES16_LFC_d[pos]))
        outfile_tiled.write('\t' + str(distance_d[pos]))
        outfile_tiled.write('\t' + str(Pfam_domain_d[pos]))
        outfile_tiled.write('\t' + str(Pfam_family_d[pos]))
        if pos in AA_cons_d:
            outfile_tiled.write('\t' + str(AA_cons_d[pos]))
        else:
            outfile_tiled.write('\t' + 'NA')
        outfile_tiled.write('\t' + str(phylo21_d[pos]))
        outfile_tiled.write('\t' + str(phast21_d[pos]))
        AAscore_list = []
        AAscore_list.append(AA_score_m2[AA_score_dm2[pos]])
        AAscore_list.append(AA_score_m1[AA_score_dm1[pos]])
        AAscore_list.append(AA_score_0[AA_score_dv2[pos]])
        AAscore_list.append(AA_score_p1[AA_score_dp1[pos]])
        AAscore_list.append(AA_score_p2[AA_score_dp2[pos]])
        AA_score = np.mean(AAscore_list)
        outfile_tiled.write('\t' + str(AA_score))
        outfile_tiled.write('\t' + str(tracrv2_rev_d[pos]))
        outfile_tiled.write('\t' + str(D16_woAA_d[pos]))
        if pos in inDelphi_infr_d:
            outfile_tiled.write('\t' + str(1-inDelphi_infr_d[pos]))
        else:
            outfile_tiled.write('\t' + 'NA')
        outfile_tiled.write('\t' + str(Bioscore_rev_d[pos]))
        outfile_tiled.write('\t' + str(VBC_score_rev_d[pos]))
        if gene not in Structure_dict:
            Structure_dict[gene] = {}
        Structure_dict[gene][pos] = RES16_LFC_d[pos]
    else:
        guide_not_found += 1
    i += 1

print(str(i) + ' guides in dataset')
print(str(guide_not_found) + ' guides skipped. missing transcript ID - no AA_cons scores')

Feature_name_list = []
Feature_name_list.append('AA_cons')
Feature_name_list.append('gRNA act')
Feature_name_list.append('In-Frames')
Feature_name_list.append('distance')
Feature_name_list.append('mod3')
Feature_name_list.append('Pfam')
Target_name = 'Abs_LFC'

data_ML = {'data':np.array(data),
                'feature_names':np.array(Feature_name_list)}

print('features')
print(data_ML['feature_names'].shape)
'''
for gene in Structure_dict:
    #os.mkdir(gene+'_Rb')
    i = 0
    x_list = []
    y_list = []
    c_list = []
    print(gene)
    for pos in sorted(Structure_dict[gene]):
        i += 1
        x_list.append(distance_d[pos])
        y_list.append(RES16_LFC_d[pos])
        c_list.append(round(Data_list_dict[pos][0], 4))
    sctr = plt.scatter(x_list, y_list, c=c_list, cmap='RdYlBu_r', vmin=0.3, linewidths=2)
    # cmap='RdYlBu_r',
    #plt.colorbar(sctr)
    plt.title(gene + 'VBC_')
    plt.xlabel('guides along cDNA sequence')
    plt.ylabel('Log 2 fold change')
    plt.savefig('VBC_D16_Okt'+'/' + gene + '_' + 'VBC_03_RB.pdf')
    plt.close()
    plt.close()
for gene in ['POLR2A']:
    #os.mkdir(gene+'_Rb')
    i = 0
    x_list = []
    y_list = []
    c_list = []
    print(gene)
    for pos in sorted(Structure_dict[gene]):
        i += 1
        x_list.append(distance_d[pos])
        y_list.append(RES16_LFC_d[pos])
        c_list.append(round(Data_list_dict[pos][0], 4))
    sctr = plt.scatter(x_list, y_list, c=c_list, cmap='RdYlBu_r', vmin=0.3, linewidths=2)
    # cmap='RdYlBu_r',
    plt.colorbar(sctr)
    plt.title(gene + 'VBC_')
    plt.xlabel('guides along cDNA sequence')
    plt.ylabel('Log 2 fold change')
    plt.savefig('VBC_D16_Okt'+'/' + gene + '_' + 'VBC_03_RB_cbar.pdf')
    plt.close()
    plt.close()

for gene in Structure_dict:
    #os.mkdir(gene+'_Rb')
    i = 0
    x_list = []
    y_list = []
    c_list = []
    print(gene)
    for pos in sorted(Structure_dict[gene]):
        i += 1
        x_list.append(distance_d[pos])
        y_list.append(RES16_LFC_d[pos])
        c_list.append(round(Data_list_dict[pos][1],4))
    sctr = plt.scatter(x_list, y_list, c=c_list, cmap='RdYlBu_r', vmin=0.2, linewidths=2)
    # cmap='RdYlBu_r',
    #plt.colorbar(sctr)
    plt.title(gene + 'Doench 2016')
    plt.xlabel('guides along cDNA sequence')
    plt.ylabel('Log 2 fold change')
    plt.savefig('VBC_D16_Okt'+'/' + gene + '_' + 'D16_RB.pdf')
    plt.close()
    plt.close()
for gene in ['POLR2A']:
    #os.mkdir(gene+'_Rb')
    i = 0
    x_list = []
    y_list = []
    c_list = []
    print(gene)
    for pos in sorted(Structure_dict[gene]):
        i += 1
        x_list.append(distance_d[pos])
        y_list.append(RES16_LFC_d[pos])
        c_list.append(round(Data_list_dict[pos][1],4))
    sctr = plt.scatter(x_list, y_list, c=c_list, cmap='RdYlBu_r', vmin=0.2, linewidths=2)
    # cmap='RdYlBu_r',
    plt.colorbar(sctr)
    plt.title(gene + 'Doench 2016')
    plt.xlabel('guides along cDNA sequence')
    plt.ylabel('Log 2 fold change')
    plt.savefig('VBC_D16_Okt'+'/' + gene + '_' + 'D16_RB_cbar.pdf')
    plt.close()
    plt.close()

'''
