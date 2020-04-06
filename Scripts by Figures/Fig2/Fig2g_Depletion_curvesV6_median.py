__author__ = 'georg.michlits'

import numpy as np
import matplotlib.pyplot as plt



##############Feb13th
#V6 - include a filter that requires all guides to have a minimum number (e.g. 1000) of wt reads in each day for the experiment
# this filter will discard 1-3 sgRNAs that are outliers in the data analysis due to increased noise level.


#
#
#
#                   TO DO!!!!
#
#
#                   or not.......
#
#
#
#



exp_type = 'n'
dict_file = open('94call_v5_mapv3_indeldict.txt','r')

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
        if guide == 'Vmn2r115_1':
            print(line)
        repl = int(column[1])
        exp_type = column[2]
        #if exp_type not in exp_type_list:
        #    print(exp_type)
        #    exp_type_list.append(exp_type)
        d = column[3].rstrip(' ')
        if 'day' in d:
            d = int(d.split('y')[1])
            mutation = column[4]
            mutation_size = mutation.split(':')[1]
            mutation_pos = mutation.split(':')[2]
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
'''
for guide in INDEL_DICT:
    for repl in INDEL_DICT[guide]:
        for day in INDEL_DICT[guide][repl][exp_type]:
            mutation = INDEL_DICT[guide][repl][exp_type]
            if mutation == 'wt:0:0':
                print(guide)
'''




print('Normalization to control guides')
Ctrls_guides = ['ROSA26_1', 'ROSA26_2', 'ROSA26_3', 'Olfr10_1', 'Olfr44_1', 'Olfr52_1']
print(Ctrls_guides)

NF = {}
for guide in Ctrls_guides:
    #print(guide)
    if guide not in NF:
        NF[guide] = {}
    for repl in sorted(INDEL_DICT[guide]):
        for day in INDEL_DICT[guide][repl]['n']:
            wt = 0
            mut = 0
            for mutation in INDEL_DICT[guide][repl]['n'][day]:
                mut_type = mutation.split(':')[0]
                count = float(INDEL_DICT[guide][repl]['n'][day][mutation])
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
#F[day] = value
F = {}
for guide in NF:
    #print(guide)
    for day in sorted(NF[guide]):
        x = np.mean(NF[guide][day])
        #print(day + '\t' + str(x))
        if day not in F:
            F[day] = []
        F[day].append(x)
print('normalization fac by days')
for day in sorted(F):
    print(str(day) + '\t' + str(1/np.mean(F[day])))
# keep option to simplify by using 3.6 average correction factor
av_fac = np.mean([8.75,2.05,1.63,2.11])
print('average correction factor ' + str(av_fac))

####################################
# FOR EACH GENE
# Determine inframe and frameshift depletion over time
# Plot count relative to wt vs time
#
#
#
####################################

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

for guide in INDEL_DICT:
    if guide not in CURVES_dict:
        CURVES_dict[guide] = {}
    x = [[], [], [], []]
    y_if = [[], [], [], []]
    y_of = [[], [], [], []]
    for repl in sorted(INDEL_DICT[guide]):
        x_index = -1
        for d in sorted(INDEL_DICT[guide][repl]['n']):
            x_index += 1
            wt = 0
            out_frame = 0
            in_frame = 0
            for mutation in INDEL_DICT[guide][repl]['n'][d]:
                mut_type = mutation.split(':')[0]
                mut_size = int(mutation.split(':')[1])
                mut_pos = int(mutation.split(':')[2])
                if -21 <= mut_pos <= 21 and mut_type in ['wt','ins','del']:
                    if mutation not in CURVES_dict[guide]:
                        CURVES_dict[guide][mutation] = {}
                        CURVES_dict[guide][mutation]['x'] = [[], [], [], []]
                        CURVES_dict[guide][mutation]['x_mean'] = []
                        CURVES_dict[guide][mutation]['y'] = [[], [], [], []]
                        CURVES_dict[guide][mutation]['y_norm'] = [[], [], [], []]
                        CURVES_dict[guide][mutation]['y_mean'] = []
                    wt = float(INDEL_DICT[guide][repl]['n'][d]['wt:0:0'])
                    count = float(INDEL_DICT[guide][repl]['n'][d][mutation])
                    if wt == 0:
                        print('ERROR - wt count = 0')
                    y = float(count/wt)*float(1/np.mean(F[d]))  #*F[d] is the normalisation factor calculated before
                    CURVES_dict[guide][mutation]['x'][x_index].append(d)
                    CURVES_dict[guide][mutation]['y'][x_index].append(y)

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
            for element in CURVES_dict[guide][mutation]['y_norm']:
                CURVES_dict[guide][mutation]['y_mean'].append(np.median(element))
            mut_type = mutation.split(':')[0]
            mut_size = int(mutation.split(':')[1])
            mut_pos = int(mutation.split(':')[2])
            if mut_type in ['del','ins']:
                x = CURVES_dict[guide][mutation]['x']
                y = CURVES_dict[guide][mutation]['y_norm']
                x_mean = CURVES_dict[guide][mutation]['x_mean']
                y_mean = CURVES_dict[guide][mutation]['y_mean']
                if mut_size > 0 and mut_size % 3 == 0:
                    plt.plot(x,y, '.b')     #I plot y not y_norm to see contibution of initial mutations
                    plt.plot(x_mean,y_mean, 'b-')
                elif mut_size == 0:
                    plt.plot(x,y, '.b')
                    plt.plot(x_mean,y_mean, 'b-')
                else:
                    plt.plot(x,y, '.r')
                    plt.plot(x_mean,y_mean, 'r-')
            plt.title(guide)
            x_labels = ['days after Cas9 induction']
    plt.tight_layout()
    plt.savefig('plots_v6-median/' + guide + '.pdf')
    plt.close()