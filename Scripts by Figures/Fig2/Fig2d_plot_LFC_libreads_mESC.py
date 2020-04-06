__author__ = 'georg.michlits'

import numpy as np
import matplotlib.pyplot as plt
import pickle

print('load lfc')
LFC_data = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Zub_n_d18_abs_LFC_d.sav','rb'))
LFC_data_2n = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Zub_2n_d18_abs_LFC_d.sav','rb'))
print('load ref count')
plasmid_count = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/mESC_C3_Zuber_d18_2n/Zub_plasmid_guide_to_count.sav','rb'))
Species = 'Ms'
CEG_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/'+Species+'/property_sav/CEG_pos_d.sav', 'rb'))
CN_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/CN_d.sav', 'rb'))

Scar_seq_guides_file = open('scarSeq_guides.txt','r')
Scar_guides = []
for line in Scar_seq_guides_file:
    guide_split = line.rstrip('\n').split('_')
    guide = guide_split[2] + '_' + guide_split[3]
    Scar_guides.append(guide)
print(Scar_guides)
print(len(Scar_guides))

print('read in Zub guides')
pos_to_REF = {}
i = 0
k = 0
ZUB_ext7_file = open('ZUB_ext7.txt','r')
for line in ZUB_ext7_file:
    if i > 0:
        column = line.rstrip('\n').split('\t')
        guidename = column[0]
        pos = column[2]
        abs_LFC = column[0]
        pos_to_REF[pos] = guidename
    i += 1
ZUB_ext7_file.close()

print(i)
print(k)

x_ls = []
y_ls = []
x_CEG = []
y_CEG = []
x_weak = []
y_weak = []
x_medium = []
y_medium = []
x_strong = []
y_strong = []

weak = []
medium = []
strong = []

LFC_screen_dict = {}
LFC_screen_dict_n = {}
sgRNA_subgroup_dict = {}

for pos in LFC_data:
    if CN_d[pos] == 1:
        guidename = pos_to_REF[pos]
        LFC_n = LFC_data[pos]
        LFC_2n = LFC_data_2n[pos]
        LFC = np.mean([LFC_n,LFC_2n])
        x_ls.append(LFC)
        y_ls.append(plasmid_count[guidename])
        if pos in CEG_d:
            x_CEG.append(LFC)
            y_CEG.append(plasmid_count[guidename])
        if guidename in Scar_guides and LFC > -3:
            LFC_screen_dict_n[guidename] = LFC_n
            LFC_screen_dict[guidename] = LFC
            sgRNA_subgroup_dict[guidename] = 'weak'
            weak.append(guidename)
            x_weak.append(LFC)
            y_weak.append(plasmid_count[guidename])
        if guidename in Scar_guides and -3 >= LFC > -5:
            LFC_screen_dict_n[guidename] = LFC_n
            LFC_screen_dict[guidename] = LFC
            sgRNA_subgroup_dict[guidename] = 'medium'
            medium.append(guidename)
            x_medium.append(LFC)
            y_medium.append(plasmid_count[guidename])
        if guidename in Scar_guides and -5 >= LFC:
            LFC_screen_dict_n[guidename] = LFC_n
            LFC_screen_dict[guidename] = LFC
            sgRNA_subgroup_dict[guidename] = 'strong'
            strong.append(guidename)
            x_strong.append(LFC)
            y_strong.append(plasmid_count[guidename])

pickle.dump(LFC_screen_dict_n,open('Scarseq_sgRNA_LFCinScreen_n_dict.sav','wb'))
pickle.dump(LFC_screen_dict,open('Scarseq_sgRNA_LFCinScreen_dict.sav','wb'))
pickle.dump(sgRNA_subgroup_dict,open('Scarseq_sgRNA_subgroups_dict.sav','wb'))

mean_CEG_depl = np.mean(x_CEG)
print(mean_CEG_depl)

print('weak' + str(len(weak)))
print(weak)
print('medium' + str(len(medium)))
print(medium)
print('strong' + str(len(strong)))
print(strong)

print(np.array(x_ls))
print(np.array(y_ls))
plt.figure(figsize=(12,10))
plt.yscale('log')
plt.scatter(x_ls,y_ls, c='grey', s=30, alpha=0.3)
plt.scatter(x_CEG,y_CEG, c='red', s=30, alpha=0.3)
plt.scatter(x_weak,y_weak, c='orange', s=50, alpha=1)
plt.scatter(x_medium,y_medium, c='mediumaquamarine', s=50, alpha=1)
plt.scatter(x_strong,y_strong, c='black', s=50, alpha=1)
plt.xlim(-12,4)
plt.tick_params(labelsize=20)
#plt.show()
plt.savefig('mESC_n18_vslibreads.png')
plt.close()
