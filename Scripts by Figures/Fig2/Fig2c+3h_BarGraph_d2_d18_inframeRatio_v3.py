__author__ = 'georg.michlits'

import numpy as np
import matplotlib.pyplot as plt
import pickle

dict_file = open('94call_v5_mapv3_indeldict.txt','r')
n_2n = 'n'
##### results come just as print out - i copy and past them into excel and prism for graph processing


print('assemble indeldict')
#indelDICT[guide][repl][type][d][mutation]=count
INDEL_DICT ={}
i = 0
dict_file.seek(0)
exp_type_list = []
for line in dict_file:
    if i > 0:
        column = line.rstrip('\n').split('\t')
        guide = column[0]
        repl = int(column[1])
        exp_type = column[2]
        #if exp_type not in exp_type_list:
        #    print(exp_type)
        #    exp_type_list.append(exp_type)
        d = column[3].rstrip(' ')
        if 'day' in d:
            d = int(d.split('y')[1])
            mutation = column[4]
            count = column[5]
            if not guide in INDEL_DICT:
                INDEL_DICT[guide] = {}
            if not repl in INDEL_DICT[guide]:
                INDEL_DICT[guide][repl] = {}
            if not exp_type in INDEL_DICT[guide][repl]:
                INDEL_DICT[guide][repl][exp_type] = {}
            if not d in INDEL_DICT[guide][repl][exp_type]:
                INDEL_DICT[guide][repl][exp_type][d] = {}
            if not mutation in INDEL_DICT[guide][repl][exp_type][d]:
                INDEL_DICT[guide][repl][exp_type][d][mutation] = count
    i += 1
print(i)

# Normalization_factor [guide] [day] = mut/wt
# we observe overall depletion of guides targeting non-essential genes (ROSA and Olfr)
# to control for that effect "cutting effect" we determine control factors

print('Normalization to control guides')
Ctrls_guides = ['ROSA26_1', 'ROSA26_2', 'ROSA26_3', 'Olfr10_1', 'Olfr44_1', 'Olfr52_1']
print(Ctrls_guides)
exp_type = n_2n
NF = {}
for guide in Ctrls_guides:
    #print(guide)
    if guide not in NF:
        NF[guide] = {}
    for repl in sorted(INDEL_DICT[guide]):
        for day in INDEL_DICT[guide][repl][n_2n]:
            wt = 0
            mut = 0
            for mutation in INDEL_DICT[guide][repl][n_2n][day]:
                mut_type = mutation.split(':')[0]
                count = float(INDEL_DICT[guide][repl][n_2n][day][mutation])
                if mut_type == 'wt':
                    wt += count
                elif mut_type == 'del':
                    mut += count
                elif mut_type == 'ins':
                    mut += count
            #now wt is count for repX dayX
            if day not in NF[guide]:
                NF[guide][day] = []
            NF[guide][day].append(float(mut/wt))
#F[day][repl] = value
F = {}
Fac = {}
for guide in NF:
    #print(guide)
    for day in sorted(NF[guide]):
        if day not in F:
            F[day] = {}
        for i in range(1,4):        # 1-4 or later i vs i-1 to convert the repl index 0,1,2 to repl value 1,2,3
            if i not in F[day]:
                F[day][i] = []          # here i is repl value
            x = NF[guide][day][i-1]     #-1 here i is repl index
            F[day][i].append(x)
print('normalization fac by days')
for day in sorted(F):
    if not day in Fac:
        Fac[day] = {}
    for repl in sorted(F[day]):
        Fac[day][repl] = np.mean(F[day][repl])
        print(str(day) + '\t' + str(Fac[day]))
# keep option to simplify by using 3.6 average correction factor
print('daily correction f average ' + str([8.75,2.05,1.63,2.11]))
av_fac = np.mean([8.75,2.05,1.63,2.11])
print('average correction factor ' + str(av_fac))


########################################################################
#     READ IN CATEGORY OF sgRNAs (Super, med, weak, ctrl)
########################################################################


sgRNA_subgroup_dict = pickle.load(open('Scarseq_sgRNA_subgroups_dict.sav','rb'))


#INDEL_DICT[guide][repl][exp_type][d][mutation] = count
# across all inframes / frameshifts
# another time individual mutations maybe

# with individual mutation. read in a dictionary that has guides and mutations
CURVES_dict = {}
'''
guide
    mutation
        x = [[],[],[],[]]
        x_mean = [[],[],[],[]]
        y_if_norm = [[],[],[],[]]
        y_of_norm = [[],[],[],[]]
        y_if_mean = [[],[],[],[]]
        y_of_mean = [[],[],[],[]]
    if mutation blabla
        add up values.
'''

Super_depleters_inframrate_d2 = []
medium_depleters_inframrate_d2 = []
weak_depleters_inframrate_d2 = []
controls_inframrate_d2 = []
Super_depleters_inframrate_d18 = []
medium_depleters_inframrate_d18 = []
weak_depleters_inframrate_d18 = []
controls_inframrate_d18 = []

print('category' + '\t' + 'guide' + '\t' + 'day2_inframes' + '\t' + 'day18_inframes')
for guide in INDEL_DICT:
    if guide not in CURVES_dict:
        CURVES_dict[guide] = {}
    x = [[],[],[],[]]
    y_if = [[],[],[],[]]
    y_of = [[],[],[],[]]
    for repl in sorted(INDEL_DICT[guide]):
        x_index = -1
        for d in sorted(INDEL_DICT[guide][repl][n_2n]):
            x_index += 1
            wt = 0
            out_frame = 0
            in_frame = 0
            for mutation in INDEL_DICT[guide][repl][n_2n][d]:
                mut_type = mutation.split(':')[0]
                mut_size = int(mutation.split(':')[1])
                mut_pos = int(mutation.split(':')[2])
                if -21 <= mut_pos <= 21 and mut_type in ['wt','ins','del']:
                    if mutation not in CURVES_dict[guide]:
                        CURVES_dict[guide][mutation] = {}
                        CURVES_dict[guide][mutation]['x'] = [[],[],[],[]]
                        CURVES_dict[guide][mutation]['x_mean'] = []
                        CURVES_dict[guide][mutation]['y'] = [[],[],[],[]]
                        CURVES_dict[guide][mutation]['y_norm'] = [[],[],[],[]]
                        CURVES_dict[guide][mutation]['y_mean'] = []
                    wt = float(INDEL_DICT[guide][repl][n_2n][d]['wt:0:0'])
                    count = float(INDEL_DICT[guide][repl][n_2n][d][mutation])
                    if wt == 0:
                        print('ERROR - wt count = 0')
                    y = float(count/wt)*float(1/Fac[d][repl])  #*F[d] is the normalisation factor calculated before
                    CURVES_dict[guide][mutation]['x'][x_index].append(d)
                    CURVES_dict[guide][mutation]['y'][x_index].append(y)

    inframes_list_d2 = []
    outframes_list_d2 = []
    inframes_list_d18 = []
    outframes_list_d18 = []
    frame_test = 'empty'
    for mutation in CURVES_dict[guide]:
        for i in range(0,4):
            day_list = [2,4,8,18]
            while len(CURVES_dict[guide][mutation]['x'][i]) < 3:
                CURVES_dict[guide][mutation]['x'][i].append(day_list[i])
                CURVES_dict[guide][mutation]['y'][i].append(0)
        y_NF = np.mean(CURVES_dict[guide][mutation]['y'][0])
        if y_NF >= 0.01:
            for i in range(0,4):
                for element in CURVES_dict[guide][mutation]['y'][i]:
                    CURVES_dict[guide][mutation]['y_norm'][i].append(element/y_NF)
            for element in CURVES_dict[guide][mutation]['x']:
                CURVES_dict[guide][mutation]['x_mean'].append(np.median(element))
            for element in CURVES_dict[guide][mutation]['y']:
                CURVES_dict[guide][mutation]['y_mean'].append(np.median(element))
            mut_type = mutation.split(':')[0]
            mut_size = int(mutation.split(':')[1])
            mut_pos = int(mutation.split(':')[2])
            if mut_type in ['del','ins']:
                x = CURVES_dict[guide][mutation]['x']
                #print(x)
                y = CURVES_dict[guide][mutation]['y']
                #print(y)
                x_mean = CURVES_dict[guide][mutation]['x_mean']
                #print(x_mean)
                y_mean = CURVES_dict[guide][mutation]['y_mean']
                #print(y_mean)
                if mut_size % 3 == 0:
                    #plt.plot(x,y, '.g')     #I plot y not y_norm to see contibution of initial mutations
                    #plt.plot(x_mean,y_mean, 'g-')
                    frame_test = 'if'
                    y_mean_d2 = y_mean[0]
                    y_mean_d18 = y_mean[3]
                else:
                    #plt.plot(x,y, '.r')
                    #plt.plot(x_mean,y_mean, 'r-')
                    frame_test = 'of'
                    y_mean_d2 = y_mean[0]
                    y_mean_d18 = y_mean[3]
        if frame_test == 'if':
            inframes_list_d2.append(y_mean_d2)
            inframes_list_d18.append(y_mean_d18)
        elif frame_test == 'of':
            outframes_list_d2.append(y_mean_d2)
            outframes_list_d18.append(y_mean_d18)
    sum_d2 = (sum(inframes_list_d2)+sum(outframes_list_d2))
    sum_d18 = (sum(inframes_list_d18)++sum(outframes_list_d18))
    if sum_d2 == 0:
        print(guide)
        sum_d2 = 0.00001
    if sum_d18 == 0:
        print(guide)
        sum_d18 = 0.00001
    guide_inf_d18 = sum(inframes_list_d18)/sum_d18
    guide_inf_d2 = sum(inframes_list_d2)/sum_d2
    if guide in sgRNA_subgroup_dict:
        if sgRNA_subgroup_dict[guide] == 'strong':
            print('strong' + '\t' + guide + '\t' + str(round(guide_inf_d2,4)) + '\t' + str(round(guide_inf_d18,4)))
            #Super_depleters_inframrate_d2.append(guide_inf_d2)
            #Super_depleters_inframrate_d18.append(guide_inf_d18)
        if sgRNA_subgroup_dict[guide] == 'medium':
            print('medium' + '\t' + guide + '\t' + str(round(guide_inf_d2,4)) + '\t' + str(round(guide_inf_d18,4)))
            #medium_depleters_inframrate_d2.append(guide_inf_d2)
            #medium_depleters_inframrate_d18.append(guide_inf_d18)
        if sgRNA_subgroup_dict[guide] == 'weak':
            print('weak' + '\t' + guide + '\t' + str(round(guide_inf_d2,4)) + '\t' + str(round(guide_inf_d18,4)))
            #weak_depleters_inframrate_d2.append(guide_inf_d2)
            #weak_depleters_inframrate_d18.append(guide_inf_d18)
    if guide in Ctrls_guides:
        print('ctrl' + '\t' + guide + '\t' + str(round(guide_inf_d2,4)) + '\t' + str(round(guide_inf_d18,4)))
        #controls_inframrate_d2.append(guide_inf_d2)
        #controls_inframrate_d18.append(guide_inf_d18)
        #plt.title(guide)
        #x_labels = ['days after Cas9 induction']
    #plt.savefig('plots_v3.1/' + guide)
    #plt.close()

