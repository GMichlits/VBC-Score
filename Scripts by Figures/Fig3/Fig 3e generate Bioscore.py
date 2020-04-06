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
Output_name = 'Bioscore_3Nov_Jun_parameters'
Feat_plot_dir = Output_name + '_Feature_plots'
Coef_plot_dir = Output_name + '_Coef_plots'
Model_dir = Output_name + '_Model'
Fit_plot_dir = Output_name + '_Fit_plots'
os.mkdir(Feat_plot_dir)
outfile_mod3pfam = open(Feat_plot_dir + '/mod3pfam_data_out.txt', 'w')
os.mkdir(Coef_plot_dir)
os.mkdir(Model_dir)
os.mkdir(Fit_plot_dir)

print('load phylo')
phylo_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/phylo21_d.sav', 'rb'))
print('load phast')
phast_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/phast21_d.sav', 'rb'))
print('load AA63')
AA_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/AA7_d.sav', 'rb'))
#print('load Doench')
#D16_woAA_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/D16_woAA_d.sav', 'rb'))
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

AA_type_cut = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/AA_type_cut_d.sav', 'rb'))
AA_type_down = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/AA_type_down_d.sav', 'rb'))
AA_type_up = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/AA_type_up_d.sav', 'rb'))

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

CEG_d = pickle.load(open('../../Gen_properties/' + Species + '/property_sav/CEG_pos_d.sav', 'rb'))
#if 'IMP' in Output_name:
#data_abs_1 = pickle.load(open('../../Gen_data/Zuber_v3/KBM7_CEG3_rel_LFC_d.sav','rb'))
#data_abs_2 = pickle.load(open('../../Gen_data/Zuber_v3/RKO_CEG3_rel_LFC_d.sav','rb'))
#data_abs_3 = pickle.load(open('../../Gen_data/Zuber_v3/MIAPACA2_CEG3_rel_LFC_d.sav','rb'))
#if 'Nov' in Output_name:
data_abs_4 = pickle.load(open('../../Gen_data/Novartis_tiled/DLD_rel90_LFC_d_m0.5.sav','rb'))
data_abs_5 = pickle.load(open('../../Gen_data/Novartis_tiled/RKO_rel90_LFC_d_m0.5.sav','rb'))
data_abs_6 = pickle.load(open('../../Gen_data/Novartis_tiled/NCIH1299_rel90_LFC_d_m0.5.sav','rb'))
#if 'Brunello' in Output_name:
#data_abs_7 = pickle.load(open('../../Gen_data/v3_DICTIONARIES/Brunello_CEG3_rel_LFC_d.sav','rb'))


i = 0
#for pos in data_abs:
#    print(pos + '\t' + str(data_abs[pos]))
#    i += 1
#    if i == 5:
#        break
all_AA = 'RHKDESTNQGAVILMFYWCP'

score_importance_train = {}
score_importance_test = {}
#IMP_3_screens = [data_abs_1, data_abs_2, data_abs_3, data_abs_4, data_abs_5, data_abs_6, data_abs_7]
IMP_3_screens = [data_abs_4, data_abs_5, data_abs_6]
x_name = ['blank','Phylo7','Phast7','Pfam','others','AA_7','AA_type','combo','random']
#x_name = ['blank','combo']
for x in x_name:
    score_importance_train[x] = []
    score_importance_test[x] = []
    i = 0
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
            if pos in phylo_d:
                if pos in AA_d and pos in AA_type_up:
                    if pos in Gene_direction_d and pos in exon_dist_start_d and pos in exon_dist_stop_d:
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
                        ##                    x Features for optimization
                        i_training += 1
                        #PHYLO
                        ['blank','Phylo7','Phast7','Pfam','others','AA_7','AA_type','combo','random']
                        if x == 'Phylo7':
                            score = phylo_d[pos]
                            data_list.append(score)
                        if x == 'Phast7':
                            score = phast_d[pos]
                            data_list.append(score)
                        if x == 'Pfam':
                            score = Pfam_domain_d[pos]
                            data_list.append(score)
                            score = Pfam_family_d[pos]
                            data_list.append(score)
                        if x == 'others':
                            score = distance_d[pos]
                            data_list.append(score)
                            data_list.append(mod3_d[pos])
                            data_list.append((exon_dist_start_d[pos]+exon_dist_stop_d[pos]))
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
                        if x == 'AA_7':
                            score = AA_d[pos]
                            data_list.append(score)
                        if x == 'AA_type':
                            data_list.append(AA_score_m2[AA_score_dm2[pos]])
                            data_list.append(AA_score_m1[AA_score_dm1[pos]])
                            data_list.append(AA_score_0[AA_score_dv2[pos]])
                            data_list.append(AA_score_p1[AA_score_dp1[pos]])
                            data_list.append(AA_score_p2[AA_score_dp2[pos]])
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

                        if x == 'combo':
                            phylo = phylo_d[pos]
                            data_list.append(phylo)
                            score = phast_d[pos]
                            data_list.append(score)
                            score = Pfam_domain_d[pos]
                            data_list.append(score)
                            score = Pfam_family_d[pos]
                            data_list.append(score)
                            score = distance_d[pos]
                            data_list.append(score)
                            data_list.append(mod3_d[pos])
                            data_list.append((exon_dist_start_d[pos]+exon_dist_stop_d[pos]))
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
                            score = AA_d[pos]
                            data_list.append(score)
                            data_list.append(AA_score_m2[AA_score_dm2[pos]])
                            data_list.append(AA_score_m1[AA_score_dm1[pos]])
                            data_list.append(AA_score_0[AA_score_dv2[pos]])
                            data_list.append(AA_score_p1[AA_score_dp1[pos]])
                            data_list.append(AA_score_p2[AA_score_dp2[pos]])
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
                        if x == 'random':
                            score = np.random.random()
                            data_list.append(score)
                        if x == 'blank':
                            score = np.random.random()
                            data_list.append(score)
                        #rand
                        #data_list.append(np.random.random())
                        # here you can insert if pos in CEG c.list append red
                        c_list.append('k')
                        #####################################################################################
                        ## y - LFC data to fit
                        data.append(data_list)
                        if -1*data_set[pos] > 0:
                            target.append(0)
                        else:
                            target.append(-1*data_set[pos])
            else:
                guide_not_found += 1
            i += 1

    print(str(i) + ' guides in dataset')
    #print(str(i_training) + ' guides in training dataset')
    #print(str(i_validation) + ' guides in validation dataset')
    #print(str(guide_found_wrong_AA_len) + ' guides found but on the edge 5pr or 3pr of CDS')
    print(str(guide_not_found) + ' guides skipped. missing transcript ID - no AA_cons scores')
    Feature_name_list = []
    ['blank','Phylo7','Phast7','Pfam','others','AA_7','AA_type','combo','random']
    if x == 'Phylo7':
        Feature_name_list.append('phylo7')
    if x == 'Phast7':
        Feature_name_list.append('Phast7')
    if x == 'AA_7':
        Feature_name_list.append('AA7')
    if x == 'Pfam':
        Feature_name_list.append('pfam_dom')
        Feature_name_list.append('pfam_fam')
    if x == 'others':
        Feature_name_list.append('dist')
        Feature_name_list.append('mod3')
        Feature_name_list.append('exlen')
        Feature_name_list.append('SA')
        Feature_name_list.append('SD')
    if x == 'AA_type':
        Feature_name_list.append('AAm2')
        Feature_name_list.append('AAm1')
        Feature_name_list.append('AA0')
        Feature_name_list.append('AAp2')
        Feature_name_list.append('AAp2')
        for letter in all_AA:
            Feature_name_list.append(letter)
    if x == 'combo':
        Feature_name_list.append('phylo7')
        Feature_name_list.append('phast7')
        Feature_name_list.append('pfam_dom')
        Feature_name_list.append('pfam_fam')
        Feature_name_list.append('dist')
        Feature_name_list.append('mod3')
        Feature_name_list.append('exlen')
        Feature_name_list.append('SA')
        Feature_name_list.append('SD')
        Feature_name_list.append('AA7')
        Feature_name_list.append('AAm2')
        Feature_name_list.append('AAm1')
        Feature_name_list.append('AA0')
        Feature_name_list.append('AAp2')
        Feature_name_list.append('AAp2')
        for letter in all_AA:
            Feature_name_list.append(letter)
    if x == 'random':
        Feature_name_list.append('random')
    if x == 'blank':
        Feature_name_list.append('random')
    #Feature_name_list.append('random2')


    Target_name = 'log2 fold change'
    c_list_feat = ['k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k','k']

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
            score_importance_train[x].append(clf.score(X_train,y_train)-blank_score_train)
            score_importance_test[x].append(clf.score(X_test,y_test)-blank_score_test)
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
    name = 'MedianFit_' + Output_name + '_' + str(i)
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
    name = 'Coef_mean_' + Output_name
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
['blank','Phylo7','Phast7','Pfam','others','AA_7','AA_type','combo','random']
c_list = ['grey','grey','grey','cadetblue','c','dodgerblue','b','grey']

score_importance_train.pop('blank') #kick out the blank value (it is substracted from all other scores)
for x in score_importance_train:
    x_name.append(x)
    y_mean.append(np.mean(score_importance_train[x]))
    y_err.append(np.std(score_importance_train[x]))
plt.bar(x_name, y_mean, color=c_list, yerr=y_err)
name = 'Train_' + Output_name
plt.savefig(name + '.pdf')
plt.close()

x_name = []
y_mean = []
y_err = []
score_importance_test.pop('blank') #kick out the blank value (it is substracted from all other scores)
for x in score_importance_test:
    x_name.append(x)
    y_mean.append(np.mean(score_importance_test[x]))
    y_err.append(np.std(score_importance_test[x]))

plt.bar(x_name, y_mean, color=c_list, yerr=y_err)
name = 'Test_' + Output_name
plt.savefig(name + '.pdf')
plt.close()
