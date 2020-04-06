
__author__ = 'georg.michlits'

import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split

Species = 'Ms'
model_input = 'Jul14_ontiled_v2.sav'
Output_name = Species + model_input.split('.')[0]
Feat_plot_dir = Output_name + '_Feature_plots'
Coef_plot_dir = Output_name + '_Coef_plots'
Model_dir = Output_name + '_Model'
Fit_plot_dir = Output_name + '_Fit_plots'
#os.mkdir(Feat_plot_dir)
#outfile_mod3pfam = open(Feat_plot_dir + '/mod3pfam_data_out.txt','w')
#os.mkdir(Coef_plot_dir)
#os.mkdir(Model_dir)
#os.mkdir(Fit_plot_dir)

print('load model')
model_VBCscore = pickle.load(open('../../Models/VBC_score/' + model_input,'rb'))
print('load Doench')
#D16_woAA_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/D16_woAA_d.sav', 'rb'))
combo_sgRNA_d = pickle.load(open('../../Gen_properties/' + Species + '/property_sav/combo_sgRNA_d.sav','rb'))
#tracrv2_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/tracrv2_3Nov_d.sav', 'rb'))
print('load inDelphi')
inDelphi_infr_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/inDelphi_infr_d.sav', 'rb'))
print('load copynumberfilter')
CN_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/CN_d.sav', 'rb'))
print('load Bioscore')
Bioscore_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/Jul14_Bioscore_d.sav', 'rb'))

###################
# NOTES
#   All dictionaries with _d contain "NAME_Property"_d = {pos : score}
#   All scores are between 0 and 1
#   pos is a guide ID string: chr12:85822107-85822137(+)
#   Dictionaries AA_cons_gn30_d has option of choosing window size "NAME_Property"_d = {pos : {1: score,
#                                                                                         3: score,
#                                                                                         5: score,.....}
#
###################
'''
Gene_direction_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/Gene_direction_d.sav', 'rb'))
exon_dist_start_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/exon_dist_start_d.sav', 'rb'))
exon_dist_stop_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/exon_dist_stop_d.sav', 'rb'))
Gene_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/Gene_d.sav', 'rb'))
UMI_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/UMI_d.sav', 'rb'))
distance_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/distance_d.sav', 'rb'))
mod3_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/mod3_d.sav', 'rb'))
Pfam_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/Pfam_all_d.sav', 'rb'))
Pfam_domain_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/Pfam_domain_d.sav', 'rb'))
Pfam_family_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/Pfam_family_d.sav', 'rb'))
Pfam_repeat_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/Pfam_repeat_d.sav', 'rb'))
ATG_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/ATG_d.sav', 'rb'))
_4T_strech_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/_4T_strech_d.sav', 'rb'))
over4T_strech_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/over4T_strech_d.sav', 'rb'))
transcript_ID_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/transcriptID_d.sav', 'rb'))
cut_frame_0_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/cut_frame_0_d.sav', 'rb'))
cut_frame_1_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/cut_frame_1_d.sav', 'rb'))
cut_frame_2_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/cut_frame_2_d.sav', 'rb'))
'''

#data_abs = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Brunello/Br_medm0.5_abs_LFC_d.sav','rb'))

infr_list = []

for pos in inDelphi_infr_d:
    infr_list.append(inDelphi_infr_d[pos])
median_infr = np.median(infr_list)
print('median inframes' + str(median_infr))

VBC_crude_score_d = {}
Bioscore_crude_score_d = {}
VBC_list = []
Bioscore_list = []
total_guides = 0
guide_excluded = 0
for pos in combo_sgRNA_d:
    data_list = []
    total_guides += 1
    #if np.random.random() <= 0.1:
    #    i_validation += 1
    #    U6_list = AA_U6PAM_d[pos]['U6']
    #    PAM_list = AA_U6PAM_d[pos]['PAM']
    #    if len(U6_list)+len(PAM_list) == 21:
    #        U6_rev = []
    #        for x in reversed(U6_list):
    #            U6_rev.append(x)
    #        list = U6_rev + PAM_list[1:]
    #        data_valid.append(list)
    #        target_valid.append(2**target_data_file[pos])
    #    else:
    #        guide_found_wrong_AA_len += 1
    #else:
    #####################################################################################
    ##                    x Features  for optimization
    data_list.append(combo_sgRNA_d[pos])
    #data_list.append(D16_woAA_d[pos])
    if pos in inDelphi_infr_d:
        data_list.append(1-inDelphi_infr_d[pos])
    else:
        data_list.append(1-median_infr)
    data_list.append(Bioscore_d[pos])
    VBC_score = round(model_VBCscore.predict([data_list])[0],5)*-1
    VBC_crude_score_d[pos] = VBC_score
    VBC_list.append(VBC_score)
    #print('VBC ' + str(VBC_score) + '\t' + 'Bioscore ' + str(Bioscore))
    if total_guides % 100000 == 0:
        print(str(total_guides/1000000) + ' million guides processed')

print('total guides: ' +str(total_guides))
print('guides_excluded: ' +str(guide_excluded))

VBC_max = np.max(VBC_list)
VBC_min = np.min(VBC_list)
#Biosc_max = np.max(Bioscore_list)
#Biosc_min = np.min(Bioscore_list)

#plt.scatter(VBC_list,Bioscore_list,alpha=0.2)
#plt.xlabel = 'VBC_score_crude'
#plt.ylabel = 'Bio_score_crude'
#plt.savefig(Output_name + '_VBCscBio_crudeScores.png')

print(VBC_max)
print(VBC_min)
#print(Biosc_max)
#print(Biosc_min)

VBC_score_d = {}
for pos in VBC_crude_score_d:
    VBC_scaled = (VBC_crude_score_d[pos]-VBC_min)/(VBC_max-VBC_min)
    #if VBC_scaled > 0.2:
    #    VBC_scaledv2 = (VBC_scaled-0.2)/0.8*0.9+0.1
    #else:
    #    VBC_scaledv2 = (VBC_scaled)/0.2*0.1
    #VBC_score_d[pos] = VBC_scaled
    VBC_score_d[pos] = VBC_scaled

#Bioscore_d = {}
#for pos in Bioscore_crude_score_d:
#    Bioscore_scaled = (Bioscore_crude_score_d[pos]-Biosc_min)/(Biosc_max-Biosc_min)
#    if Bioscore_scaled > 0.6:
#        Bio_scaledv2 = (Bioscore_scaled-0.6)/0.4*0.9+0.1
#    else:
#        Bio_scaledv2 = (Bioscore_scaled)/0.6*0.1
#    Bioscore_d[pos] = Bio_scaledv2

#pickle.dump(VBC_score_d, open(Species + model_input.split('.')[0] +'_d.sav','wb'))
#pickle.dump(Bioscore_d, open('Bioscore_3IMPsc_CEG3_rel_v2_d.sav','wb'))