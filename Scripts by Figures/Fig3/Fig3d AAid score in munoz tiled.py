__author__ = 'georg.michlits'

import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split

outfile_data = open('AAid_Nov20_outfile_data.txt','w')

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

AA_type_cut = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/AA_type_cut_d.sav', 'rb'))
AA_type_down = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/AA_type_down_d.sav', 'rb'))
AA_type_up = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/AA_type_up_d.sav', 'rb'))

print('load other useful dictionaries')
Gene_d = pickle.load(open('../../Gen_properties/' + Species + '/property_sav/Gene_d.sav', 'rb'))
nt30_d = pickle.load(open('../../Gen_properties/' + Species + '/property_sav/nt30_d.sav', 'rb'))
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

folder = 'AAid_Nov20/'
#polar = ['R','H','K','D','E','S','T','N','Q']
#apolar = ['G','A','V','I','L','M','F','Y','W']
#pecial = ['C','P']

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

all_AA = 'RHKDESTNQGAVILMFYWCP'
#AA_letters = 'RHKDESTNQGAVILMFYWCP'
AA_letters = '00000000000000000000'
coeff_13window = {
'R':-0.023366,
'H':-0.018707,
'K':0.015808,
'D':0.003499,
'E':-0.002382,
'S':-0.013784,
'T':-0.027023,
'N':0.010490,
'Q':-0.004988,
'G':-0.020781,
'A':-0.014429,
'V':-0.029615,
'I':-0.031432,
'L':-0.057271,
'M':-0.004342,
'F':-0.009628,
'Y':-0.042177,
'W':-0.040596,
'C':-0.038575,
'P':0.008981}

outfile_data.write('gene' + '\t' + 'sgRNA pos' + '\t' + 'nt30seq' + '\t' + 'LFC - DLD1' + '\t' + "5' to 3' position" + '\t' + 'AA id score')

for gene in Structure_dict:
    #os.mkdir(gene+'_Rb')
    x_list = []
    y_list = []
    c_list = []
    for pos in sorted(Structure_dict[gene]):
        i += 1
        LFC = DLD1_LFC_d[pos]
        dist = distance_d[pos]
        x_list.append(dist)
        y_list.append(LFC)
        gene = Gene_d[pos]
        nt30 = nt30_d[pos]
        data_list = []
        data_list.append(AA_score_m2[AA_score_dm2[pos]])
        data_list.append(AA_score_m1[AA_score_dm1[pos]])
        data_list.append(AA_score_0[AA_score_dv2[pos]])
        data_list.append(AA_score_p1[AA_score_dp1[pos]])
        data_list.append(AA_score_p2[AA_score_dp2[pos]])
        AA_score = (data_list[0]*-0.003297 + data_list[1]*0.000640 + data_list[2] * -0.004437 + data_list[3] * -0.002845 + data_list[4] * -0.003812)
        binary_13_w_score = 0
        if pos in AA_type_cut:
            for letter in all_AA:
                UP_down_range = 6 # 6 corres[onds to window 13: cutsite 1 + 6 up  +6 down.
                if letter in AA_type_cut[pos]:
                    binary_13_w_score += coeff_13window[letter]
                elif letter in AA_type_up[pos][0:UP_down_range]:
                    binary_13_w_score += coeff_13window[letter]
                elif letter in AA_type_down[pos][0:UP_down_range]:
                    binary_13_w_score += coeff_13window[letter]
                else:
                    binary_13_w_score += 0
        else:
            for letter in AA_letters:
                binary_13_w_score = 0
        AA_identity_score = (AA_score + binary_13_w_score) * -1
        outfile_data.write('\n' + gene + '\t' + pos + '\t' + nt30 + '\t' + str(LFC) + '\t' + str(dist) + '\t' + str(AA_identity_score))
        c_list.append(AA_identity_score)
    name = 'DLD1' + gene
    sctr = plt.scatter(x_list, y_list, c=c_list, cmap='RdYlBu_r', linewidth=2)
    plt.colorbar(sctr)
    plt.xlabel('guides along coding DNA sequence')
    plt.ylabel('enrichment LFC')
    plt.savefig(folder + 'V2AAScore' + name + '.pdf')
    plt.close()
    plt.close()
    plt.close()