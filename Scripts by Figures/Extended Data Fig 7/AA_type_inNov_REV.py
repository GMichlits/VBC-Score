__author__ = 'georg.michlits'

import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split


print('load dataset')
DLD1_LFC_d = pickle.load(open('../../Gen_data/Novartis_tiled/DLD_abs_LFC_d.sav', 'rb'))

Species = 'Hm'
print('load Features')
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

print('load other useful dictionaries')
Gene_d = pickle.load(open('../../Gen_properties/' + Species + '/property_sav/Gene_d.sav', 'rb'))
distance_d = pickle.load(open('../../Gen_properties/' + Species + '/property_sav/distance_d.sav', 'rb'))
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

print('data')
i = 0
guide_not_found = 0
data = []
for pos in DLD1_LFC_d:
    if pos in AA_score_dp1 and pos in Gene_d:
        gene = Gene_d[pos]
        if not gene in Structure_dict:
            Structure_dict[gene] = {}
        Structure_dict[gene][pos] = DLD1_LFC_d[pos]
    else:
        guide_not_found += 1
    i += 1

print(str(i) + ' guides in dataset')
print(str(guide_not_found) + ' guides skipped. missing transcript ID - no AA_cons scores')

folder = 'AA_groups_v2REV/'
polar = ['R','H','K','D','E','S','T','N','Q']
apolar = ['G','A','V','I','L','M','F','Y','W']
special = ['C','P']

# for gene in Structure_dict:
#     #os.mkdir(gene+'_Rb')
#     x_list = []
#     y_list = []
#     c_list = []
#     for pos in sorted(Structure_dict[gene]):
#         i += 1
#         x_list.append(i)
#         y_list.append(DLD1_LFC_d[pos])
#         if AA_score_dv2[pos] in polar:
#             c_list.append('purple')
#         if AA_score_dv2[pos] in apolar:
#             c_list.append('g')
#         if AA_score_dv2[pos] in special:
#             c_list.append('y')
#
#     name = 'DLD1' + gene
#     sctr = plt.scatter(x_list, y_list, c=c_list)
#     plt.xlabel(name)
#     plt.ylabel('enrichment LFC')
#     plt.savefig(folder +'group' + name + '.pdf')
#     print(name + 'group.pdf')
#     plt.close()
#     plt.close()
#     plt.close()

for gene in Structure_dict:
    #os.mkdir(gene+'_Rb')
    x_list = []
    y_list = []
    c_list = []
    for pos in sorted(Structure_dict[gene]):
        i += 1
        x_list.append(distance_d[pos])
        y_list.append(DLD1_LFC_d[pos])
        data_list = []
        data_list.append(AA_score_m2[AA_score_dm2[pos]])
        data_list.append(AA_score_m1[AA_score_dm1[pos]])
        data_list.append(AA_score_0[AA_score_dv2[pos]])
        data_list.append(AA_score_p1[AA_score_dp1[pos]])
        data_list.append(AA_score_p2[AA_score_dp2[pos]])
        AA_score = np.mean(data_list)
        c_list.append(AA_score)
    name = 'DLD1' + gene
    sctr = plt.scatter(x_list, y_list, c=c_list, cmap='RdYlBu_r', linewidth=2)
    plt.colorbar(sctr)
    plt.xlabel('guides along coding DNA sequence')
    plt.ylabel('enrichment LFC')
    plt.savefig(folder + 'V2AAScore' + name + '.pdf')
    plt.close()
    plt.close()
    plt.close()