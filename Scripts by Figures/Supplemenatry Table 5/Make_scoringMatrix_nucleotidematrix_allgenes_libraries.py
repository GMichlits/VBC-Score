__author__ = 'georg.michlits'

import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split



#model_name = '3IMPsc_medm1_rel_v2'
#### add all improved tracr version screens

Screen_data_namelist = ['../../Gen_data/v3_DICTIONARIES/TKOv3_HAP1d18_rel_LFC_d.sav']
Species = 'Hm'
Screen_data_namelist.append('../../Gen_data/Yusa/Exp_data/Exp0abs_LFC_d.sav')
Screen_data_namelist.append('../../Gen_data/Yusa/Exp_data/Exp1abs_LFC_d.sav')
Screen_data_namelist.append('../../Gen_data/Yusa/Exp_data/Exp12abs_LFC_d.sav')
Screen_data_namelist.append('../../Gen_data/v3_DICTIONARIES/Brunello_rel_LFC_d.sav')
Screen_data_namelist.append('../../Gen_data/v3_DICTIONARIES/TKOv1_HTC116_rel_LFC_d.sav')
Screen_data_namelist.append('../../Gen_data/v3_DICTIONARIES/Wang15_Raji_Exp2rel_LFC_d.sav')
Screen_data_namelist.append('../../Gen_data/v3_DICTIONARIES/Wang17_MOLM13_Exp136rel_LFC_d.sav')
Screen_data_namelist.append('../../Gen_data/v3_DICTIONARIES/GeckoV2_Exp_12_HT29_LARGE_INTESTINE_rel_LFC_d.sav')
Screen_data_namelist.append('../../Gen_data/v3_DICTIONARIES/Avana_HCC1143_genCR_Exp227rel_LFC_d.sav')


Screen_data_namelist.append('../../Gen_data/Zuber_v3/MIAPACA2_rel_LFC_d.sav')
'''
Screen_data_namelist = []
Species = 'Ms'
Screen_data_namelist.append('../../Gen_data/v3_DICTIONARIES/CrUMI_2n_rel_LFC_d.sav')
Screen_data_namelist.append('../../Gen_data/v3_DICTIONARIES/Zub_2n_d18_rel_LFC_d.sav')
'''

#print(model_name + ' ' + data_name)

##### HartMoffat TKOv1
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/TKOv1_HTC116_CEG3_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/TKOv1_HTC116_abs_LFC_d.sav', 'rb'))
##### HartMoffat TKOv3
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/TKOv3_HAP1d18_CEG3_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/TKOv3_HAP1d18_abs_d.sav', 'rb'))
##### Wang15_Raji
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Wang15_Raji_Exp2CEG3_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Wang15_Raji_Exp2abs_LFC_d.sav', 'rb'))
##### Wang15_K562
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Wang15_K562_Exp3CEG3_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Wang15_K562_Exp3abs_LFC_d.sav', 'rb'))
##### Wang17_MOLM13
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Wang17_MOLM13_Exp136CEG3_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Wang17_MOLM13_Exp136abs_LFC_d.sav', 'rb'))
#####Wang17_PL21
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Wang17_PL21_Exp128CEG3_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Wang17_PL21_Exp128abs_LFC_d.sav', 'rb'))
##### Geckov2_HT29
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/GeckoV2_Exp_12_HT29_LARGE_INTESTINE_CEG3_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/GeckoV2_Exp_12_HT29_LARGE_INTESTINE_abs_d.sav', 'rb'))
##### Geckov2_K562
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/GeckoV2_Exp_13_K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE_CEG3_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/GeckoV2_Exp_13_K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE_abs_d.sav', 'rb'))
##### Geckov2 PC3
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Gecko_v2_Exp_25_PC3_PROSTATE_CEG3_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Gecko_v2_Exp_25_PC3_PROSTATE_abs_d.sav', 'rb'))

##### Yusa HL60 Exp3
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Yusa/Exp_data_thirdCEG/Exp3third_abs_LFC.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Yusa/Exp_data/Exp3abs_LFC_d.sav', 'rb'))
##### Yusa HT1080 Exp4
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Yusa/Exp_data_thirdCEG/Exp4third_abs_LFC.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Yusa/Exp_data/Exp4abs_LFC_d.sav', 'rb'))
##### Avana Karpas
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Avana_Karpas_genCR_Exp285_CEG3_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Avana_Karpas_genCR_Exp285_abs_LFC_d.sav', 'rb'))
##### Avana HCC1143 - Exp83
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Avana_HCC1143_genCR_Exp227CEG3_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Avana_HCC1143_genCR_Exp227abs_LFC_d.sav', 'rb'))
#### Brunello DATA LFC-0.5
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Brunello/Br_medm0.5_abs_LFC_d.sav', 'rb'))
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Brunello/Br_medm1_abs_LFC_d.sav', 'rb'))
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Brunello/Br_medm1.5_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Brunello/Brunello_abs_LFC_d.sav', 'rb'))
#### BRunello CEG3
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Brunello_CEG3_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Brunello_abs_LFC_d.sav', 'rb'))
#### IMP KBM7 DATA
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/KBM7_medm1_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/KBM7_abs_LFC_d.sav', 'rb'))
##### IMP RKO
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/RKO_medm1_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/RKO_abs_LFC_d.sav', 'rb'))
##### IMP MIApaca2
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/MIAPACA2_medm1_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/MIAPACA2_abs_LFC_d.sav', 'rb'))
##### CrUMI mESC Data_2n
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/CrUMI_2n_CEG3_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/CrUMI_2n_abs_LFC_d.sav', 'rb'))
##### CrUMI mESC Data_n
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/CrUMI_n_CEG3_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/CrUMI_n_abs_LFC_d.sav', 'rb'))
##### IMP mESC Data_2n
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Zub_2n_d18_CEG3_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Zub_2n_d18_abs_LFC_d.sav', 'rb'))
##### IMP mESC Data_n
#Screen_data_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Zub_n_d18_CEG3_abs_LFC_d.sav', 'rb'))
#Screen_data_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Zub_n_d18_abs_LFC_d.sav', 'rb'))



#print('load model_name')
#model_T1_AAw5_6prop_abs_LFC = pickle.load(open('/Users/georg.michlits/Desktop/Models/VBC_score/model_Nov_rel.sav','rb'))
# gn30 Ms is gn63 Hm
#print('load Doench_with AApositional information')
#D16_AA_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/D16_AA_d.sav', 'rb'))
#print('load Doench_without AApositional information')
#D16_woAA_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/D16_woAA_d.sav', 'rb'))
#print('load inDelphi')
#inDelphi = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/inDelphi_infr_d.sav', 'rb'))

#inDelphi_fr_d = {}
#for pos in inDelphi:
#    inDelphi_fr_d[pos] = 1-inDelphi[pos]

#VBC_score_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/VBC_score_3IMPsc_CEG3_rel_v2_d.sav', 'rb'))
#VBC_score_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/VBC_score_KBM7_medm1_rel_v2_d.sav', 'rb'))
#VBC_score_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/VBC_score_KBM7_CEG3_rel_v2_d.sav', 'rb'))
#VBC_score_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/VBC_score_Nov_rel_v2_d.sav', 'rb'))
#VBC_score_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/VBC_score_3Nov_rel_v2_d.sav', 'rb'))
#VBC_score_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/VBC_score_3IMPsc_medm1_abs_v2_d.sav', 'rb'))
#VBC_score_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/VBC_score_3IMPsc_CEG3_abs_v2_d.sav', 'rb'))
#VBC_score_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/VBC_score_KBM7_medm1_abs_v2_d.sav', 'rb'))
#VBC_score_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/VBC_score_KBM7_CEG3_abs_v2_d.sav', 'rb'))



###################  Load basic information required
#print('load CEG')
#CEG_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/' + Species + '/property_sav/CEG_pos_d.sav', 'rb'))
print('load Gene_d')
Gene_d = pickle.load(open('../../Gen_properties/' + Species + '/property_sav/Gene_d.sav', 'rb'))
print('load nt_30')
nt30_d = pickle.load(open('../../Gen_properties/' + Species + '/property_sav/nt30_d.sav', 'rb'))

###########################################################################

def extract_bestworst_pos(Screen_data_rel_d, Gene_d):
    Struc_dict = {}
    best_worst = {}
    for pos in Screen_data_rel_d:
        if pos in Gene_d:
            rel_LFC = Screen_data_rel_d[pos]
            gene = Gene_d[pos]
            if not gene in Struc_dict:
                Struc_dict[gene] = {'pos': [], 'rel_LFC': []}
            Struc_dict[gene]['pos'].append(pos)
            Struc_dict[gene]['rel_LFC'].append(rel_LFC)
    for gene in Struc_dict:
        gene_pos_list = Struc_dict[gene]['pos']
        index_max = Struc_dict[gene]['rel_LFC'].index(np.max(Struc_dict[gene]['rel_LFC']))
        index_min = Struc_dict[gene]['rel_LFC'].index(np.min(Struc_dict[gene]['rel_LFC']))
        best_pos = gene_pos_list[index_max]
        worst_pos = gene_pos_list[index_min]
        best_worst[gene] = {'b':best_pos, 'w':worst_pos}
    return (best_worst)

score_matrix_output = open(Species + '_scorematrix_out_checkLibrcomp.txt','w')
kcluster_matrix_output = open(Species + '_into_kcluster_checkLibrcomp.txt','w')
kcluster_matrix_output.write('Screen_Name\tdata')

for Screen_data_name in Screen_data_namelist:
    print(Screen_data_name + 'load')
    Screen_data_rel = pickle.load(open(Screen_data_name,'rb'))
    print(Screen_data_name + 'processed')

    score_matrix_input_bias = {}
    for i in range(0, 30):
        score_matrix_input_bias[i] = {}
        score_matrix_input_bias[i]['A'] = 0
        score_matrix_input_bias[i]['C'] = 0
        score_matrix_input_bias[i]['G'] = 0
        score_matrix_input_bias[i]['T'] = 0

    for pos in Screen_data_rel:
        nt30_seq = nt30_d[pos]
        for i in range(0, 30):
            letter = nt30_seq[i]
            score_matrix_input_bias[i][letter] += 1

    score_matrix_raw = {}
    for i in range(0, 30):
        score_matrix_raw[i] = {}
        score_matrix_raw[i]['A'] = 0
        score_matrix_raw[i]['C'] = 0
        score_matrix_raw[i]['G'] = 0
        score_matrix_raw[i]['T'] = 0

    best_worst = extract_bestworst_pos(Screen_data_rel, Gene_d)
    for gene in best_worst:
        best_pos = best_worst[gene]['b']
        worst_pos = best_worst[gene]['w']
        best_nt30_seq = nt30_d[best_pos]
        worst_nt30_seq = nt30_d[worst_pos]
        for i in range(0, 30):
            best_letter = best_nt30_seq[i]
            worst_letter = worst_nt30_seq[i]
            score_matrix_raw[i][best_letter] += 1
            score_matrix_raw[i][worst_letter] += -1

    score_matrix_output.write('\n' + Screen_data_name)
    kcluster_matrix_output.write('\n' + Screen_data_name)
    score_matrix_output.write('\n' + '\t' + 'A' + '\t' + 'C' + '\t' + 'G' + '\t' + 'T' + '\t' + 'A' + '\t' + 'C' + '\t' + 'G' + '\t' + 'T')
    for i in range(0,30):
        score_matrix_output.write('\n' + str(i+1) + '\t' + str(score_matrix_raw[i]['A']))
        score_matrix_output.write('\t' + str(score_matrix_raw[i]['C']))
        score_matrix_output.write('\t' + str(score_matrix_raw[i]['G']))
        score_matrix_output.write('\t' + str(score_matrix_raw[i]['T']))
        score_matrix_output.write('\t' + str(score_matrix_input_bias[i]['A']))
        score_matrix_output.write('\t' + str(score_matrix_input_bias[i]['C']))
        score_matrix_output.write('\t' + str(score_matrix_input_bias[i]['G']))
        score_matrix_output.write('\t' + str(score_matrix_input_bias[i]['T']))
        score_current_position = {}
        for letter in ['A','C','G','T']:
            if score_matrix_input_bias[i][letter] == 0:
                score_current_position[letter] = 0
            else:
                score_current_position[letter] = score_matrix_raw[i][letter]/score_matrix_input_bias[i][letter]
            kcluster_matrix_output.write('\t' + str(score_current_position[letter]))
