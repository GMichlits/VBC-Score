__author__ = 'georg.michlits'

import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

#scaler = StandardScaler()
#scaler.fit(X_train)
#X_train = scaler.transform(X_train)
#X_test = scaler.transform(X_test)

Species = 'Hm'
Output_name = 'singleAA_3tiled_multiINs_LARGEwindows_v2'
Feat_plot_dir = Output_name + '_Feature_plots'
Coef_plot_dir = Output_name + '_Coef_plots'
Model_dir = Output_name + '_Model'
Fit_plot_dir = Output_name + '_Fit_plots'
os.mkdir(Feat_plot_dir)
outfile_mod3pfam = open(Feat_plot_dir + '/mod3pfam_data_out.txt', 'w')
os.mkdir(Coef_plot_dir)
os.mkdir(Model_dir)
os.mkdir(Fit_plot_dir)

#print('load phylo')
#phylo_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/phylo_v2_d.sav', 'rb'))
#print('load phast')
#phast_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/phast_v2_d.sav', 'rb'))
#print('load AA63')
#AA_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_cons_gn63_d.sav', 'rb'))
#print('load Doench')

#AA_score_dv2 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_residue_dv2.sav', 'rb'))
#AA_score_dp1 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_residue_dp1.sav', 'rb'))
#AA_score_dp2 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_residue_dp2.sav', 'rb'))
#AA_score_dm1 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_residue_dm1.sav', 'rb'))
#AA_score_dm2 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_residue_dm2.sav', 'rb'))

#AA_resIN1_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_resIN1_v3_d.sav', 'rb'))
#AA_resIN3_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_resIN3_v3_d.sav', 'rb'))
#AA_resIN5_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_resIN5_v3_d.sav', 'rb'))
#AA_resIN7_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_resIN7_v3_d.sav', 'rb'))
#AA_resIN9_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_resIN9_v3_d.sav', 'rb'))
#AA_resIN11_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_resIN11_v3_d.sav', 'rb'))

AA_type_cut = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_type_cut_d.sav','rb'))
AA_type_down = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_type_down_d.sav','rb'))
AA_type_up = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_type_up_d.sav','rb'))

#phylo_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/phylo21_d.sav', 'rb'))
#AA_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA7_d.sav', 'rb'))
#inDelphi_infr_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/inDelphi_infr_d.sav', 'rb'))
#distance_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/distance_d.sav', 'rb'))
#mod3_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/mod3_d.sav', 'rb'))
#Pfam_domain_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Pfam_domain_d.sav', 'rb'))
#Pfam_family_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Pfam_family_d.sav', 'rb'))
#exon_dist_start_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/exon_dist_start_d.sav', 'rb'))
#exon_dist_stop_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/exon_dist_stop_d.sav', 'rb'))
#D16_woAA_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/D16_woAA_d.sav', 'rb'))
#tracrv2_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/tracrv2_mms_d.sav', 'rb'))
#tracrv2_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/tracrv2_mhm_d.sav', 'rb'))
#tracrv2_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/tracrv2_d.sav', 'rb'))
#IMPv2_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/IMP22_d.sav', 'rb'))
#print('load inDelphi')
#inDelphi_infr_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/inDelphi_infr_d.sav', 'rb'))
#print('load dist')
#distance_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/distance_d.sav', 'rb'))
#mod3_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/mod3_d.sav', 'rb'))
#print('load pfam')
#Pfam_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Pfam_all_d.sav', 'rb'))
#Pfam_domain_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Pfam_domain_d.sav', 'rb'))
#Pfam_family_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Pfam_family_d.sav', 'rb'))
#Pfam_repeat_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Pfam_repeat_d.sav', 'rb'))
#print('load T_streches')
#_4T_strech_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/_4T_strech_d.sav', 'rb'))
#over4T_strech_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/over4T_strech_d.sav', 'rb'))
#print('load cutting_frame')
#cut_frame_0_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/cut_frame_0_d.sav', 'rb'))
#cut_frame_1_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/cut_frame_1_d.sav', 'rb'))
#cut_frame_2_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/cut_frame_2_d.sav', 'rb'))
#print('load splice')
#Gene_direction_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/Gene_direction_d.sav', 'rb'))
#exon_dist_start_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/exon_dist_start_d.sav', 'rb'))
#exon_dist_stop_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/exon_dist_stop_d.sav', 'rb'))

#print('load gaps')
#AA_gaps_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/AA_gaps_d.sav', 'rb'))

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


AA_scorev2 = {'H' : 3.40143101,
              'R' : 1.08052375,
            'K' : -4.453883147,
            'D' : -7.705176369,
            'E' : -2.049637515,
            'S' : -5.176921052,
            'T' : 3.145321586,
            'N' : -5.88211047,
            'Q' : -2.337120609,
            'G' : 3.364400709,
            'A' : -5.766952958,
            'V' : 3.855088089,
            'I' : 3.994762418,
            'L' : 2.750942054,
            'M' : 1.06067543,
            'F' : -1.795480355,
            'Y' : 8.083789415,
            'W' : 13.8994159,
            'C' : 6.413725606,
            'P' : -1.248312289,}


AA_score = {'R':0.529061666,
            'H':3.945168241,
            'K':-2.898681179,
            'D':-6.521588809,
            'E':1.007236132,
            'S':-4.564723986,
            'T':2.82254932,
            'N':-3.11570683,
            'Q':-2.831809165,
            'G':1.738132035,
            'A':-4.948095132,
            'V':4.282030392,
            'I':4.970418996,
            'L':3.096945294,
            'M':-0.339088887,
            'F':2.707394875,
            'Y':8.343544002,
            'W':9.588584727,
            'C':-1.213702742,
            'P':-4.574936124,}

CEG_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/' + Species + '/property_sav/CEG_pos_d.sav', 'rb'))
#data_abs_1 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/KBM7_medm1_rel_LFC_d.sav','rb'))
#data_abs_2 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/RKO_medm1_rel_LFC_d.sav','rb'))
#data_abs_3 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Zuber_v3/MIAPACA2_medm1_rel_LFC_d.sav','rb'))

data_abs_1 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Novartis_tiled/DLD_rel90_LFC_d_m0.5.sav','rb'))
data_abs_2 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Novartis_tiled/RKO_rel90_LFC_d_m0.5.sav','rb'))
data_abs_3 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Novartis_tiled/NCIH1299_rel90_LFC_d_m0.5.sav','rb'))

#data_abs_1 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Zub_2n_d18_CEG3_rel_LFC_d.sav', 'rb'))
#data_abs_2 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Zub_n_d18_CEG3_rel_LFC_d.sav', 'rb'))

#data_abs_7 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Brunello_CEG3_rel_LFC_d.sav', 'rb'))
#data_abs_8 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Yusa/Exp_data_thirdCEG/Exp3third_rel_LFC_d.sav', 'rb'))
#data_abs_9 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Yusa/Exp_data_thirdCEG/Exp4third_rel_LFC_d.sav', 'rb'))
#data_abs_10 = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/TKOv3_HAP1d18_CEG3_rel_LFC_d.sav', 'rb'))

#data_abs = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/Novartis_tiled/rel90_LFC_d_m0.5.sav','rb'))
i = 0
#for pos in data_abs:
#    print(pos + '\t' + str(data_abs[pos]))
#    i += 1
#    if i == 5:
#        break

unpolar = ['W','Y','L','I','V','F']
all_AA = 'RHKDESTNQAVILMFYWCGP'

score_importance_train = {}
score_importance_test = {}
x_name = ['1','3','5','7','9','11','13','15','17','19','21','23','29']
IMP_3_screens = [data_abs_1, data_abs_2, data_abs_3]#,data_abs_4, data_abs_5, data_abs_6,data_abs_7, data_abs_8, data_abs_9,data_abs_10]
#Nov3_screens = [data_abs_1n, data_abs_2n, data_abs_3n]#,data_abs_4, data_abs_5, data_abs_6,data_abs_7, data_abs_8, data_abs_9,data_abs_10]
for x in x_name:
    score_importance_train[x] = []
    score_importance_test[x] = []
    i = 0
    i_n = 0
    i_training = 0
    #i_validation = 0
    guide_not_found = 0
    guide_found_wrong_AA_len = 0
    data = []
    c_list = []
    #data_valid = []
    target = []
    mod3_0 = []
    mod3_1 = []
    Pfam_0 = []
    Pfam_1 = []
    #target_valid = []
    for data_set in IMP_3_screens:
        for pos in data_set:
            if pos in AA_type_cut and pos in AA_type_down and pos in AA_type_up:
                data_list = []
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
                i_training += 1
                #PHYLO
                [1,3,5,7,9,11,13,15,17,19,21,23,29]
                #if x == 'tracrv2':
                    #data_list.append(tracrv2_d[pos])
                #if x == 'D16':
                    #data_list.append(D16_woAA_d[pos])
                #if x == 'AA_hit':
                 #   for letter in all_AA:
                 #       if letter == AA_score_d[pos]:
                  #          data_list.append(1)
                  #      else:
                   #         data_list.append(0)
                    #data_list.append(AA_score[AA_score_d[pos]])
                for letter in all_AA:
                    UP_down_range = int((int(x)-1)/2)
                    if letter in AA_type_cut[pos]:
                        data_list.append(1)
                    elif letter in AA_type_up[pos][0:UP_down_range]:
                        data_list.append(1)
                    elif letter in AA_type_down[pos][0:UP_down_range]:
                        data_list.append(1)
                    else:
                        data_list.append(0)
                #if x == 'comb':
                    #data_list.append(IMPv2_d[pos])
                    #data_list.append(tracrv2_d[pos])
                    #data_list.append(D16_woAA_d[pos])
                    #data_list.append(AA_score[AA_score_d[pos]])

                #if x == 'inDelphi':
                    #data_list.append(1-inDelphi_infr_d[pos])
                #if x == 'sgRNAscore+inDelphi':
                    #data_list.append(tracrv2_d[pos])
                    #data_list.append(D16_woAA_d[pos])
                    #data_list.append(AA_score_d[pos])
                    #data_list.append(1-inDelphi_infr_d[pos])
                if x == '':
                    score = np.random.random()
                    data_list.append(score)

                #rand
                #data_list.append(np.random.random())
                # here you can insert if pos in CEG c.list append red
                c_list.append('k')
                #####################################################################################
                ## y - LFC data to fit
                data.append(data_list)
                if (-1*data_set[pos]) > 0:
                    target.append(0)
                else:
                    target.append(-1*data_set[pos])
            i += 1

    print(str(i) + ' guides in dataset IMP')
    print(str(i_n) + ' guides in dataset Nov')
    #print(str(i_training) + ' guides in training dataset')
    #print(str(i_validation) + ' guides in validation dataset')
    #print(str(guide_found_wrong_AA_len) + ' guides found but on the edge 5pr or 3pr of CDS')
    print(str(guide_not_found) + ' guides skipped. missing transcript ID - no AA_cons scores')
    Feature_name_list = []
    ['tracrv2','D16','AA_hit','A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','res_in3','rand']
    if x == 'IMPv2':
        Feature_name_list.append('IMPv2')
    if x == 'tracrv2':
        Feature_name_list.append('tracrv2')
    if x == 'D16':
        Feature_name_list.append('Rule Set 2')
    if x in ['1','3','5','7','9','11','13','15','17','19','21','23','29']:
        for letter in all_AA:
            Feature_name_list.append(letter)
        #Feature_name_list.append('AA_hit')
    for letter in all_AA:
        if x == letter:
            Feature_name_list.append(letter)
    if x == 'res_IN1':
        for letter in all_AA:
            Feature_name_list.append(letter)
    if x == 'res_IN3':
        for letter in all_AA:
            Feature_name_list.append(letter)
    if x == 'res_IN5':
        for letter in all_AA:
            Feature_name_list.append(letter)
    if x == 'res_IN7':
        for letter in all_AA:
            Feature_name_list.append(letter)
    if x == 'res_IN9':
        for letter in all_AA:
            Feature_name_list.append(letter)
    if x == 'res_IN11':
        for letter in all_AA:
            Feature_name_list.append(letter)
    if x == 'comb':
        #Feature_name_list.append('IMPv2')
        Feature_name_list.append('tracrv2')
        Feature_name_list.append('Rule Set 2')
        Feature_name_list.append('Hart2017')
    if x == 'sgRNAscore-Zub':
        Feature_name_list.append('tracrv2')
        Feature_name_list.append('Rule Set 2')
        Feature_name_list.append('Hart2017')
    if x == 'inDelphi':
        Feature_name_list.append('inDelphi')
    if x == 'sgRNAscore+inDelphi':
        Feature_name_list.append('IMPv2')
        Feature_name_list.append('tracrv2')
        Feature_name_list.append('Rule Set 2')
        Feature_name_list.append('Hart2017')
        Feature_name_list.append('inDelphi')
    if x == '':
        Feature_name_list.append('random')
    #Feature_name_list.append('random2')


    Target_name = 'log2 fold change'
    c_list_feat = ['k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k']

    data_ML = {'data': np.array(data),
                'feature_names': np.array(Feature_name_list),
                'target_names': np.array(Target_name),
                'target': np.array(target)}

    print('Data input')
    print(data_ML['data'].shape)
    print('features')
    print(data_ML['feature_names'].shape)
    print('target')
    print(data_ML['target'].shape)

    i = 0
    for index, feature_name in enumerate(data_ML['feature_names']):
        plt.title(feature_name)
        plt.scatter(data_ML['target'], data_ML['data'][:,index], c=c_list, s=10, alpha=0.4)
        gradient, intercept, r_value, p_value, std_err = stats.linregress(data_ML['target'], data_ML['data'][:,index])
        plt.xlabel(Target_name, size=15)
        plt.ylabel(feature_name, size=15)
        mn = np.min(data_ML['target'])
        mx = np.max(data_ML['target'])
        x1 = np.linspace(mn, mx, 500)
        y1 = gradient*x1+intercept
        plt.plot(x1, y1, c_list_feat[i], linewidth=2)
        plt.plot(x1, y1, linewidth=2)
        i += 1
        plt.savefig(Feat_plot_dir + '/' + feature_name + '.png')
        plt.close()
        plt.close()
        plt.close()

    ###make mod3 plot
    outfile_mod3pfam.write('mod3_0\tmod3_1\tpfam_0\tpfam_1')
    n = max([len(mod3_0), len(mod3_1), len(Pfam_0), len(Pfam_1)])
    for i in range(n):
        if i < len(mod3_0):
            outfile_mod3pfam.write('\n' + str(mod3_0[i]))
        else:
            outfile_mod3pfam.write('\n')
        if i < len(mod3_1):
            outfile_mod3pfam.write('\t' + str(mod3_1[i]))
        else:
            outfile_mod3pfam.write('\t')
        if i < len(Pfam_0):
            outfile_mod3pfam.write('\t' + str(Pfam_0[i]))
        else:
            outfile_mod3pfam.write('\t')
        if i < len(Pfam_1):
            outfile_mod3pfam.write('\t' + str(Pfam_1[i]))
        else:
            outfile_mod3pfam.write('\t')

    RMS_list = []
    clf_score_list = []
    Test_score_list = []
    coef_list = []
    n_runs = 9
    print('current analysis: ' + x)
    for i in range(n_runs):
        X_train, X_test, y_train, y_test = train_test_split(data_ML['data'], data_ML['target'])

        #scaler = StandardScaler()
        #scaler.fit(X_train)

        #X_train = scaler.transform(X_train)
        #X_test = scaler.transform(X_test)

        clf = LinearRegression()
        clf.fit(X_train, y_train)

        print('clf_score train-test ' + str(clf.score(X_train,y_train)) + '-' + str(clf.score(X_test,y_test)))
        if not x == 'blank':
            print(x)
            score_importance_train[x].append(clf.score(X_train,y_train))
            score_importance_test[x].append(clf.score(X_test,y_test))
        else:
            score_importance_train[x].append(clf.score(X_train,y_train))
            score_importance_test[x].append(clf.score(X_test,y_test))

        clf_score_list.append('clf_score train-test ' + str(clf.score(X_train,y_train)) + '-' + str(clf.score(X_test,y_test)))

        predicted = clf.predict(X_test)
        expected = y_test
        plt.scatter(expected ,-1*predicted, alpha=0.4, s=10, c=c_list) # ULIS favourite style plt.scatter((np.log2(expected)), (1-(predicted)))
        plt.axis('tight')
        #plt.ylim(-4,-2)
        #plt.xlim(-10,0)
        #plt.plot([-10,0], [-10,0], '--k')
        plt.xlabel('True LFC-meassure', size=15)
        plt.ylabel('New score', size=15)
        gradient, intercept, r_value, p_value, std_err = stats.linregress(expected,-1*predicted)
        mn = np.min(expected)
        mx = np.max(expected)
        x1 = np.linspace(mn,mx,500)
        y1 = gradient*x1+intercept
        plt.plot(x1,y1, c='r', linewidth=2)
        plt.tight_layout()
        name = 'Fit_' + Output_name + '_' + str(i)
        plt.savefig(Fit_plot_dir + '/' + name + '.png')
        #plt.xlim(0,1)
        plt.close()

        print("RMS: %s" % np.sqrt(np.mean((predicted - expected) ** 2)))
        RMS_list.append(np.sqrt(np.mean((predicted - expected) ** 2)))

        Test_score_list.append(clf.score(X_test,y_test))
        coef_list.append(clf.coef_)
        plt.bar(Feature_name_list,clf.coef_)
        name = 'Coef_' + Output_name + '_' + str(i)
        plt.savefig(Coef_plot_dir + '/' + name + '.pdf')
        plt.close()
        model_name = 'model_' + Output_name + '_' + str(i)
        scaler_name = 'scaler_' + Output_name + '_' + str(i)
        pickle.dump(clf, open(Model_dir + '/' +model_name + '.sav', 'wb'))
        #pickle.dump(scaler, open(Model_dir + '/' +scaler_name + '.sav', 'wb'))

    i = Test_score_list.index(np.median(Test_score_list))
    print('Median performance model:')
    print('run_' + str(i))
    print(clf_score_list[i])
    print('load_model_' + str(i))
    median_model = pickle.load(open(Model_dir + '/model_' + Output_name + '_' + str(i) + '.sav', 'rb'))
    mm_predicted = median_model.predict(data_ML['data'])
    mm_expected = data_ML['target']
    median_model.score(data_ML['data'], data_ML['target'])


    plt.scatter(mm_expected, -1*mm_predicted, alpha=0.4, s=10, c=c_list) # ULIS favourite style plt.scatter((np.log2(expected)), (1-(predicted)))
    plt.axis('tight')
    #plt.ylim(-4,-2)
    #plt.xlim(-10,0)
    #plt.plot([-10,0], [-10,0], '--k')
    plt.xlabel('log2 fold change', size=15)
    plt.ylabel('New score', size=15)
    plt.tight_layout()
    name = 'MedianFit_' + Output_name + x + '_' + str(i)
    gradient, intercept, r_value, p_value, std_err = stats.linregress(mm_expected,-1*mm_predicted)
    mn = np.min(mm_expected)
    mx = np.max(mm_expected)
    x1 = np.linspace(mn,mx,500)
    y1 = gradient*x1+intercept
    plt.plot(x1,y1, c='r', linewidth=2)
    plt.savefig(name + '.png')
    #plt.xlim(0,1)
    plt.close()

    coef_reshape = []
    for element in Feature_name_list:
        coef_reshape.append([])

    for i in range(len(Feature_name_list)):
        for element in range(n_runs):
            coef_reshape[i].append(-1*coef_list[element][i])

    coef_mean_list = []
    coef_std_list = []
    for element in coef_reshape:
        coef_mean_list.append(np.mean(element))
        coef_std_list.append(np.std(element))
    plt.bar(Feature_name_list, coef_mean_list, color=c_list_feat, yerr=coef_std_list)
    name = 'Coef_mean_' + Output_name + x
    plt.savefig(name + '.png')
    plt.close()
    print(coef_mean_list)

    if x == 'blank':
        print('generate blank scores')
        blank_score_train = np.mean(score_importance_train['blank'])
        blank_score_test = np.mean(score_importance_test['blank'])
        print('blank_train: ' + str(blank_score_train) + '\t' + 'blank_test: ' + str(blank_score_test))


x_name = []
y_mean = []
y_err = []
['res_IN1','res_IN3','res_IN5','res_IN7','res_IN9','res_IN11','','b','b5']
c_list = ['deepskyblue','deepskyblue','deepskyblue','deepskyblue','deepskyblue','deepskyblue','deepskyblue','deepskyblue','deepskyblue','deepskyblue','deepskyblue','deepskyblue','deepskyblue','tomato','orangered','orangered','orangered','orangered','r','r','r','r','r','r','saddlebrown','saddlebrown','saddlebrown','brown','brown','grey','k','k','k']
variable = Output_name
#score_importance_train.pop('blank') #kick out the blank value (it is substracted from all other scores)
print('printout_individual clf.scores on train_data_set')
for x in score_importance_train:
    x_name.append(x)
    y_mean.append(np.mean(score_importance_train[x]))
    y_err.append(np.std(score_importance_train[x]))
    print(x + '\t' + str(y_mean) + '\t' + str(y_err))
plt.figure(figsize=(4,2))
plt.bar(x_name, y_mean, color=c_list, yerr=y_err)
name = variable + '_train_' + Output_name
plt.tick_params(labelsize=16)
#plt.ylim(0,0.12)
plt.tight_layout()
plt.savefig(name + '.pdf')
plt.close()


x_name = []
y_mean = []
y_err = []
#score_importance_test.pop('blank') #kick out the blank value (it is substracted from all other scores)
print('printout_individual clf.scores on test_data_set')
for x in score_importance_test:
    x_name.append(x)
    y_mean.append(np.mean(score_importance_test[x]))
    y_err.append(np.std(score_importance_test[x]))
    print(x + '\t' + str(y_mean) + '\t' + str(y_err))

plt.figure(figsize=(4,2))
plt.bar(x_name, y_mean, color=c_list, yerr=y_err)
plt.tick_params(labelsize=16)
#plt.ylim(0,0.12)
plt.tight_layout()
name = variable + '_test_' + Output_name
plt.savefig(name + '.pdf')
plt.close()
