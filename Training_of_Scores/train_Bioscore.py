
__author__ = 'georg.michlits'

import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split

Species = 'Ms'
model_input = '3Nov.sav'
Output_name = 'INVESTIGATE1'
Feat_plot_dir = Output_name + '_Feature_plots'
Coef_plot_dir = Output_name + '_Coef_plots'
Model_dir = Output_name + '_Model'
Fit_plot_dir = Output_name + '_Fit_plots'
os.mkdir(Feat_plot_dir)
outfile_mod3pfam = open(Feat_plot_dir + '/mod3pfam_data_out.txt','w')
os.mkdir(Coef_plot_dir)
os.mkdir(Model_dir)
os.mkdir(Fit_plot_dir)

print('load model')
model_Bioscore = pickle.load(open('../../Models/Bioscore/Junmodel_' + model_input,'rb'))
print('load phylo')
phylo_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/phylo21_d.sav', 'rb'))
print('load phast')
phast_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/phast21_d.sav', 'rb'))
print('load AA63')
if Species == 'Ms':
    AA_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/AA7_d.sav', 'rb'))
else:
    AA_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/AA7_d.sav', 'rb'))
print('load Doench')
D16_woAA_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/D16_woAA_d.sav', 'rb'))
#print('load inDelphi')
#inDelphi_infr_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/inDelphi_infr_d.sav', 'rb'))
print('load dist')
distance_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/distance_d.sav', 'rb'))
mod3_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/mod3_d.sav', 'rb'))
print('load pfam')
#Pfam_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Pfam_all_d.sav', 'rb'))
Pfam_domain_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/Pfam_domain_d.sav', 'rb'))
Pfam_family_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/Pfam_family_d.sav', 'rb'))
#Pfam_repeat_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Pfam_repeat_d.sav', 'rb'))
print('load copynumberfilter')
CN_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/CN_d.sav', 'rb'))
#print('load T_streches')
#_4T_strech_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/_4T_strech_d.sav', 'rb'))
#over4T_strech_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/over4T_strech_d.sav', 'rb'))
#print('load cutting_frame')
#cut_frame_0_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/cut_frame_0_d.sav', 'rb'))
#cut_frame_1_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/cut_frame_1_d.sav', 'rb'))
#cut_frame_2_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/cut_frame_2_d.sav', 'rb'))
print('load splice')
Gene_direction_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/Gene_direction_d.sav', 'rb'))
exon_dist_start_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/exon_dist_start_d.sav', 'rb'))
exon_dist_stop_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/exon_dist_stop_d.sav', 'rb'))
print('load AAtype')
AA_score_dm2 = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/AA_residue_d-2.sav', 'rb'))
AA_score_dm1 = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/AA_residue_d-1.sav', 'rb'))
AA_score_dv2 = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/AA_residue_dv2.sav', 'rb'))
AA_score_dp1 = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/AA_residue_d+1.sav', 'rb'))
AA_score_dp2 = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/AA_residue_d+2.sav', 'rb'))
print('load AAtype scoretables')
AA_score_m2 = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/AA_score_m2.sav', 'rb'))
AA_score_m1 = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/AA_score_m1.sav', 'rb'))
AA_score_0 = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/AA_score_0.sav', 'rb'))
AA_score_p1 = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/AA_score_p1.sav', 'rb'))
AA_score_p2 = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/AA_score_p2.sav', 'rb'))
print('load AA for 20 features in AA_window')
AA_type_cut = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/AA_type_cut_d.sav', 'rb'))
AA_type_down = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/AA_type_down_d.sav', 'rb'))
AA_type_up = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/AA_type_up_d.sav', 'rb'))

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


all_AA = 'RHKDESTNQGAVILMFYWCP'

phylo_list = []
phast_list = []
AA_list = []
Pfam_family_blank = 0.5
Pfam_domain_blank = 0.5
distance_blank = 0.5
mod3_blank = 0.5
exlen_list = []
SA_blank = 21
SD_blank = 22
AA_score = 0
#AA_letters = 'RHKDESTNQGAVILMFYWCP'
AA_letters = '00000000000000000000'

for pos in phylo_d:
    if not np.isnan(phylo_d[pos]):
        phylo_list.append(phylo_d[pos])
median_phylo_list = np.median(phylo_list)
print('median phylo ' + str(median_phylo_list))
for pos in phast_d:
    if not np.isnan(phast_d[pos]):
        phast_list.append(phast_d[pos])
median_phast_list = np.median(phast_list)
print('median phast ' + str(median_phast_list))
for pos in AA_d:
    if not np.isnan(AA_d[pos]):
        AA_list.append(AA_d[pos])
if len(AA_d) > 0:
    median_AA_list = np.median(AA_list)
else:
    median_AA_list = 0
print('median AA ' + str(median_AA_list))
for pos in exon_dist_stop_d:
    exlen_list.append(exon_dist_stop_d[pos]+exon_dist_start_d[pos])
median_exlen_list = np.median(exlen_list)
print('median exlen ' + str(median_exlen_list))


VBC_crude_score_d = {}
Bioscore_crude_score_d = {}
VBC_list = []
Bioscore_list = []
total_guides = 0
guide_excluded = 0
for pos in D16_woAA_d:
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
    if pos in phylo_d:
        if not np.isnan(phylo_d[pos]):
            data_list.append(phylo_d[pos])
        else:
            data_list.append(median_phylo_list)
    else:
        data_list.append(median_phylo_list)
    if pos in phast_d:
        if not np.isnan(phast_d[pos]):
            data_list.append(phast_d[pos])
        else:
            data_list.append(median_phast_list)
    else:
        data_list.append(median_phast_list)
    if pos in Pfam_domain_d:
        data_list.append(Pfam_domain_d[pos])
    else:
        data_list.append(0.5)
    if pos in Pfam_family_d:
        data_list.append(Pfam_family_d[pos])
    else:
        data_list.append(0.5)
    if pos in distance_d:
        data_list.append(distance_d[pos])
    else:
        data_list.append(0.5)
    if pos in mod3_d:
        data_list.append(mod3_d[pos])
    else:
        data_list.append(0.5)
    if pos in exon_dist_start_d:
        data_list.append((exon_dist_start_d[pos]+exon_dist_stop_d[pos]))
    else:
        data_list.append(median_exlen_list)
    if pos in Gene_direction_d and pos in exon_dist_start_d:
        if Gene_direction_d[pos] == '+':
            score = exon_dist_start_d[pos]
            if score > 20:
                data_list.append(21)
            else:
                data_list.append(score)
        else:
            score = exon_dist_stop_d[pos]
            if score > 20:
                data_list.append(21)
            else:
                data_list.append(score)
        if Gene_direction_d[pos] == '-':
            score = exon_dist_start_d[pos]
            if score > 20:
                data_list.append(21)
            else:
                data_list.append(score)
        else:
            score = exon_dist_stop_d[pos]
            if score > 20:
                data_list.append(21)
            else:
                data_list.append(score)
    else:
        data_list.append(21)
        data_list.append(21)
    if pos in AA_d:
        if not np.isnan(AA_d[pos]):
            if 0 <= AA_d[pos] <= 1:
                data_list.append(AA_d[pos])
            else:
                data_list.append(median_AA_list)
        else:
            data_list.append(median_AA_list)
    else:
        data_list.append(median_AA_list)
    if pos in AA_score_dv2:
        data_list.append(AA_score_m2[AA_score_dm2[pos]])
        data_list.append(AA_score_m1[AA_score_dm1[pos]])
        data_list.append(AA_score_0[AA_score_dv2[pos]])
        data_list.append(AA_score_p1[AA_score_dp1[pos]])
        data_list.append(AA_score_p2[AA_score_dp2[pos]])
    else:
        data_list.append(0)
        data_list.append(0)
        data_list.append(0)
        data_list.append(0)
        data_list.append(0)
    if pos in AA_type_cut:
        for letter in all_AA:
            UP_down_range = 6 # 6 corres[onds to window 13: cutsite 1 + 6 up  +6 down.
            if letter in AA_type_cut[pos]:
                data_list.append(1)
            elif letter in AA_type_up[pos][0:UP_down_range]:
                data_list.append(1)
            elif letter in AA_type_down[pos][0:UP_down_range]:
                data_list.append(1)
            else:
                data_list.append(0)
    else:
        for letter in AA_letters:
            data_list.append(0)
    for element in data_list:
        if np.isnan(element):
            test = 'nan'
            print('nan found')
    Bioscore = round(model_Bioscore.predict([data_list])[0],5)*-1
    Bioscore_crude_score_d[pos] = Bioscore
    Bioscore_list.append(Bioscore)
    #score = tracrv2_d[pos]
    #data_list.append(score)
    #score = 1-inDelphi_infr_d[pos]
    #data_list.append(score)
    #test = 'ok'
    #for element in data_list:
    #    if np.isnan(element):
    #        test = 'nan'
    #if not 1 >= AA_d[pos][7] >= 0:
    #    test = 'nan'
    #if test == 'ok':
    #    VBC_score = round(model_Bioscore.predict([data_list])[0],5)*-1
    #    VBC_crude_score_d[pos] = VBC_score
    #    VBC_list.append(VBC_score)
    #print('VBC ' + str(VBC_score) + '\t' + 'Bioscore ' + str(Bioscore))
    if total_guides % 100000 == 0:
        print(str(total_guides/1000000) + ' million guides processed')

print('total guides: ' + str(total_guides))
print('guides_excluded: ' + str(guide_excluded))

#VBC_max = np.max(VBC_list)
#VBC_min = np.min(VBC_list)
Biosc_max = np.max(Bioscore_list)
Biosc_min = np.min(Bioscore_list)

#plt.scatter(VBC_list,Bioscore_list,alpha=0.2)
#plt.xlabel = 'VBC_score_crude'
#plt.ylabel = 'Bio_score_crude'
#plt.savefig(Output_name + '_VBCscBio_crudeScores.png')

#print(VBC_max)
#print(VBC_min)
print(Biosc_max)
print(Biosc_min)

#VBC_score_d = {}
#for pos in VBC_crude_score_d:
#    VBC_scaled = (VBC_crude_score_d[pos]-VBC_min)/(VBC_max-VBC_min)
#    if VBC_scaled > 0.6:
#        VBC_scaledv2 = (VBC_scaled-0.6)/0.4*0.9+0.1
#    else:
#        VBC_scaledv2 = (VBC_scaled)/0.6*0.1
#    VBC_score_d[pos] = VBC_scaledv2

Bioscore_d = {}
Bioscore_scaled_list = []
for pos in Bioscore_crude_score_d:
    Bioscore_scaled = (Bioscore_crude_score_d[pos]-Biosc_min)/(Biosc_max-Biosc_min)
    if Bioscore_scaled > 0.5:
        Bio_scaledv2 = (Bioscore_scaled-0.5)/0.5*0.9+0.1
        Bioscore_scaled_list.append(Bio_scaledv2)
    else:
        Bio_scaledv2 = (Bioscore_scaled)/0.5*0.1
        Bioscore_scaled_list.append(Bio_scaledv2)
    Bioscore_d[pos] = Bio_scaledv2

#pickle.dump(VBC_score_d, open('HM_XXX_VBC_score_' + model_input.split('.')[0] +'_v2_d.sav','wb'))
pickle.dump(Bioscore_d, open(Species + '_Jul14_' + model_input.split('.')[0] + '_3Nov_d.sav','wb'))