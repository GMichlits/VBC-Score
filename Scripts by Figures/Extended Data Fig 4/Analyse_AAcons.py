__author__ = 'georg.michlits'

import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split


print('load dataset')
RES16_LFC_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Novartis_tiled/DLD_abs_LFC_d.sav', 'rb'))

Species = 'Hm'
print('load Features')
D16_AA_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/D16_AA_d.sav', 'rb'))
AA_cons_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA7_d.sav', 'rb'))
distance_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/distance_d.sav', 'rb'))
mod3_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/mod3_d.sav', 'rb'))
Pfam_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Pfam_all_d.sav', 'rb'))
inDelphi_infr_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/inDelphi_infr_d.sav', 'rb'))
AA_score_dp1 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_residue_d+1.sav', 'rb'))
print('load other useful dictionaries')
Gene_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/' + Species + '/property_sav/Gene_d.sav', 'rb'))
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

print('generatefeatures for model')
i = 0
guide_not_found = 0
data = []
for pos in RES16_LFC_d:
    if pos in AA_cons_d:
        if pos in inDelphi_infr_d:
            data_list = []
            #####################################################################################
            ##                    x Features for MODEL
            AA_cons = AA_cons_d[pos]
            data_list.append(AA_cons)
            D16_AA = D16_AA_d[pos]
            data_list.append(D16_AA)
            inDelphi_infr = inDelphi_infr_d[pos]
            data_list.append(1-inDelphi_infr)
            Distance = distance_d[pos]
            data_list.append(Distance)
            mod3 = mod3_d[pos]
            data_list.append(mod3)
            Pfam = Pfam_all_d[pos]
            data_list.append(Pfam)

            gene = Gene_d[pos]
            Data_list_dict[pos] = data_list
            if gene not in Structure_dict:
                Structure_dict[gene] = {}
            Structure_dict[gene][pos] = RES16_LFC_d[pos]
        else:
            guide_not_found += 1
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

for gene in Structure_dict:
    #os.mkdir(gene+'_Rb')
    i = 0
    x_list = []
    y_list = []
    c_list = []
    for pos in sorted(Structure_dict[gene]):
        i += 1
        x_list.append(distance_d[pos])
        y_list.append(RES16_LFC_d[pos])
        c_list.append(round(Data_list_dict[pos][0],4))

    sctr = plt.scatter(x_list, y_list, c=c_list, cmap='RdYlBu_r', linewidths=2)
    plt.colorbar(sctr)
    plt.xlabel('guides along coding DNA sequence')
    plt.ylabel('enrichment LFC')
    plt.savefig('AA_cons'+'/' + gene + '_'+ 'RB_r.pdf')
    plt.close()
    plt.close()