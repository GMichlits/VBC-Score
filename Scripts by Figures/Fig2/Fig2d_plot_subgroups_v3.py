
__author__ = 'georg.michlits'

guidedepl_file = open('guide_depl_MD.txt','r')
import numpy as np
import pickle
import matplotlib.pyplot as plt

mut_2 = pickle.load(open('mutations_dict_2_top30_relALL_n.sav','rb'))
mut_18 = pickle.load(open('mutations_dict_18_top30_relALL_n.sav','rb'))

name = 'all'

w = 22
weak = ['Dcaf13_2', 'Mrpl33_5', 'Rhox2f_1', 'Ints10_2', 'Vmn2r115_1', 'Chek1_1', 'Traip_3', 'Clns1a_1', 'Ints2_1', 'Taf8_5', 'Eif3b_5', 'Eif3b_4', 'Doxl2_4', 'Slc30a9_4', 'Jak1_2', 'Mrpl33_6', 'Doxl2_6', 'Nup93_5', 'Vmn2r115_6', 'Slc30a9_2', 'Traip_4', 'Rhox2f_3']
m = 11
medium = ['Psmd2_3', 'Stat3_4', 'Sub1_1', 'Snapc5_1', 'Nup214_4', 'Hscb_3', 'Ncbp1_2', 'Ints2_5', 'Eif3a_5', 'Nup214_2', 'Rack1_1']
s = 42
strong = ['Txn1_1', 'Socs3_4', 'Ruvbl2_1', 'Rad51_4', 'Stat3_2', 'Kif11_5', 'Jak1_5', 'Psma2_3', 'Kif18b_2', 'Socs3_2', 'Spc25_2', 'Rad51_1', 'Traip_1', 'Psma1_4', 'Nelfb_6', 'Ncbp1_5', 'Rrm1_2', 'Hscb_1', 'Traip_6', 'Taf8_4', 'Rack1_2', 'Kif18b_4', 'Clns1a_5', 'Cdk1_4', 'Psma2_2', 'Ggps1_3', 'Chek1_5', 'Psmd2_4', 'Ruvbl2_2', 'Eif3a_4', 'Srf_3', 'Rrm1_4', 'Cdk1_2', 'Dcaf13_4', 'Psma1_3', 'Rpl13a_1', 'Nelfb_2', 'Ggps1_2', 'Spc25_6', 'Snapc5_3', 'Txn1_2', 'Rpl13a_4']

#weak = ['Srf_2','Rhox2f_1','Traip_2','Eif3b_5','Hscb_3','Dcaf13_2','Nup214_4','Traip_3','Clns1a_1','Jak1_2','Traip_4','Taf8_5','Sub1_2','Chek1_1']
#medium = ['Jak1_5','Snapc5_3','Eif3a_5','Cdk1_2','Rhox2f_3','Taf8_4','Cdk1_2','Sub1_1','Dcaf13_4','Clns1a_5','Rpl35a_4','Rpl35a_6','Snapc5_1','Nup214_2','Eif3b_4']
#strong = ['Ruvbl2_1','Rpl13a_1','Psmd2_4','Nelfb_6','Rpl13a_4','Kif18b_2','Rad51_4','Psma1_4','Kif11_5','Kif18b_4','Ggps1_3','Eif3a_4','Rad51_1','Cdk1_4','Psma2_2','Nelfb_2','Spc25_6','Socs3_2','Traip_1','Rrm1_4','Ggps1_2','Socs3_4','Srf_3','Kif11_1','Ruvbl2_2','Ncbp1_5','Psma2_3','Txn1_2','Rrm1_2','Txn1_1','Stat3_4','Traip_6','Stat3_2','Rack1_1','Spc25_2','Psma1_3','Ncbp1_2','Rack1_2','Chek1_5','Psmd2_3','Hscb_1','Rbx1_1']
ctrl = ['ROSA26_1','ROSA26_2','ROSA26_3','Olfr10_1','Olfr52_1','Olfr44_1']
all = weak + medium + strong
merge_guides_indel = {}
Subgroup = {}
Subgroup['ctrl'] = ctrl
Subgroup['weak'] = weak
Subgroup['medium'] = medium
Subgroup['strong'] = strong
Subgroup['all'] = all


for guide in mut_2:
    if guide in Subgroup[name]:
        for mutation in mut_2[guide]:
            type = mutation.split(':')[0]
            dist = mutation.split(':')[1]
            pos_start = mutation.split(':')[2]
            if type == 'del':
                if -21 < int(pos_start) < 21:
                    x = -1*int(dist)
                    #print(mutation + '\t' + str(mut_2[guide][mutation]))
                    if '1' in mut_2[guide][mutation]:
                        repl_1_value = mut_2[guide][mutation]['1']
                    else:
                        repl_1_value = 0.0
                    if '2' in mut_2[guide][mutation]:
                        repl_2_value = mut_2[guide][mutation]['2']
                    else:
                        repl_2_value = 0.0
                    if repl_2_value > 5:
                        print(guide + '\t' + mutation + '\t' + str(mut_2[guide][mutation]))
                    if '3' in mut_2[guide][mutation]:
                        repl_3_value = mut_2[guide][mutation]['3']
                    else:
                        repl_3_value = 0.0
                    y = np.array([repl_1_value,repl_2_value,repl_3_value]) #3 replicates
                    if x not in merge_guides_indel:
                        merge_guides_indel[x] = np.array([0.0,0.0,0.0])
                    merge_guides_indel[x] += y
            if type == 'ins':
                if -21 < int(pos_start) < 21:
                    x = 1*int(dist)
                    #print(mutation + '\t' + str(mut_2[guide][mutation]))
                    if '1' in mut_2[guide][mutation]:
                        repl_1_value = mut_2[guide][mutation]['1']
                    else:
                        repl_1_value = 0.0
                    if '2' in mut_2[guide][mutation]:
                        repl_2_value = mut_2[guide][mutation]['2']
                    else:
                        repl_2_value = 0.0
                    if repl_2_value > 5:
                        print(guide + '\t' + mutation + '\t' + str(mut_2[guide][mutation]))
                    if '3' in mut_2[guide][mutation]:
                        repl_3_value = mut_2[guide][mutation]['3']
                    else:
                        repl_3_value = 0.0
                    y = np.array([repl_1_value,repl_2_value,repl_3_value]) #3 replicates
                    if x not in merge_guides_indel:
                        merge_guides_indel[x] = np.array([0.0,0.0,0.0])
                    merge_guides_indel[x] += y
merge_guides_del_day2 = merge_guides_indel

for x in sorted(merge_guides_indel):
    print('day2' + '\t' + str(x) + '\t' + str(merge_guides_indel[x]))


merge_guides_indel = {}
for guide in mut_18:
    if guide in Subgroup[name]:
        for mutation in mut_18[guide]:
            type = mutation.split(':')[0]
            dist = mutation.split(':')[1]
            pos_start = mutation.split(':')[2]
            if type == 'del':
                if -21 < int(pos_start) < 21:
                    x = -1*int(dist)
                    #print(mutation + '\t' + str(mut_2[guide][mutation]))
                    if '1' in mut_18[guide][mutation]:
                        repl_1_value = mut_18[guide][mutation]['1']
                    else:
                        repl_1_value = 0.0
                    if '2' in mut_18[guide][mutation]:
                        repl_2_value = mut_18[guide][mutation]['2']
                    else:
                        repl_2_value = 0.0
                    if repl_2_value > 5:
                        print(guide + '\t' + mutation + '\t' + str(mut_18[guide][mutation]))
                    if '3' in mut_18[guide][mutation]:
                        repl_3_value = mut_18[guide][mutation]['3']
                    else:
                        repl_3_value = 0.0
                    y = np.array([repl_1_value,repl_2_value,repl_3_value]) #3 replicates
                    if x not in merge_guides_indel:
                        merge_guides_indel[x] = np.array([0.0,0.0,0.0])
                    merge_guides_indel[x] += y
            if type == 'ins':
                if -21 < int(pos_start) < 21:
                    x = 1*int(dist)
                    #if dist == '1':
                    #    print(guide + '\t' + mutation + '\t' +str(mut_18[guide][mutation]))
                    #print(mutation + '\t' + str(mut_2[guide][mutation]))
                    if '1' in mut_18[guide][mutation]:
                        repl_1_value = mut_18[guide][mutation]['1']
                    else:
                        repl_1_value = 0.0
                    if '2' in mut_18[guide][mutation]:
                        repl_2_value = mut_18[guide][mutation]['2']
                    else:
                        repl_2_value = 0.0
                    if repl_2_value > 5:
                        print(guide + '\t' + mutation + '\t' + str(mut_18[guide][mutation]))
                    if '3' in mut_18[guide][mutation]:
                        repl_3_value = mut_18[guide][mutation]['3']
                    else:
                        repl_3_value = 0.0
                    y = np.array([repl_1_value,repl_2_value,repl_3_value]) #3 replicates
                    if x not in merge_guides_indel:
                        merge_guides_indel[x] = np.array([0.0,0.0,0.0])
                    merge_guides_indel[x] += y
merge_guides_del_day18 = merge_guides_indel

for x in sorted(merge_guides_indel):
    print('day18' + '\t' + str(x) + '\t' + str(merge_guides_indel[x]))

merge_guides_del_day18[0] = np.array([0.0,0.0,0.0])
merge_guides_del_day2[0] = np.array([0.0,0.0,0.0])

x_list = []
y_d2 = []
y_d2err = []
y_d18 = []
y_d18err = []
c_list = []
for x in range(-21,7):
    x_list.append(x)
    if x in merge_guides_del_day2:
        y_d2.append(np.mean(merge_guides_del_day2[x]))
        y_d2err.append(np.std(merge_guides_del_day2[x]))
    else:
        y_d2.append(0.0)
        y_d2err.append(0.0)
    if x in merge_guides_del_day18:
        y_d18.append(-1*np.mean(merge_guides_del_day18[x]))
        y_d18err.append(np.std(merge_guides_del_day18[x]))
    else:
        y_d18.append(0.0)
        y_d18err.append(0.0)
    if x % 3 == 0:
        c_list.append('b')
    else:
        c_list.append('r')

## normalize bar size so that total bar area = 1  - brackat ''' away if not wanted
'''
sumyd2 = np.sum(y_d2)
sumyd18 = -1*np.sum(y_d18)
for index in range(len(y_d2)):
    y_d2[index] = y_d2[index]/sumyd2
    y_d18[index] = y_d18[index]/sumyd18
    y_d2err[index] = y_d2err[index]/sumyd2
    y_d18err[index] = y_d18err[index]/sumyd18
'''
##################

N = len(x_list)
width = 0.75       # the width of the bars

fig, ax = plt.subplots(1,1,figsize=(14,7))
        #rects1 = ax.bar(ind, men_means, width, color='r', yerr=men_std)
ylim = max(max(y_d2),min(y_d18)*-1)*1.1
plt.ylim(ylim*-1,ylim)
plt.axhline(0, color = 'black')
plt.axvline(0, color = 'black')
rects_positive_day2 = ax.bar(x_list, y_d2, width, color = c_list, yerr=y_d2err)
rects_negative_day18 = ax.bar(x_list, y_d18, width, color = c_list, yerr=y_d18err)
        #women_means = (25, 32, 34, 20, 25)
        #women_std = (3, 5, 2, 3, 3)
        #rects2 = ax.bar(ind + width, women_means, width, color='y', yerr=women_std)
        # add some text for labels, title and axes ticks

#ax.set_title(type + 'day2vs18_range+-' +str(pos_range_value) + 'bp' )
#ax.set_ylabel(type + '_freq')
ax.tick_params(labelsize=20)
#ax.set_xlabel(type + '_size [bp]')
plt.tight_layout()
#plt.show()
plt.savefig('v5.3_data/' + name + '_v3_d2vsd18_top30_relALL_n.pdf')
plt.close()