__author__ = 'georg.michlits'

import numpy as np
import matplotlib.pyplot as plt
import pickle

print('load lfc screen data')
LFC_data = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Zub_n_d18_abs_LFC_d.sav','rb'))
LFC_data_2n = pickle.load(open('/Users/georg.michlits/Desktop/Gen_data/v3_DICTIONARIES/Zub_2n_d18_abs_LFC_d.sav','rb'))
print('load copynumbers')
CN_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/CN_d.sav', 'rb'))
print('load pfam')
Pfam_all_d = pickle.load(open('/Users/georg.michlits/Desktop/Gen_properties/Ms/property_sav/Pfam_all_d.sav', 'rb'))


pfam_groups = {}
pfam_groups['neutral'] = {'no_pfam':[], 'yes_pfam':[]}
pfam_groups['weak'] = {'no_pfam':[], 'yes_pfam':[]}
pfam_groups['medium'] = {'no_pfam':[], 'yes_pfam':[]}
pfam_groups['strong'] = {'no_pfam':[], 'yes_pfam':[]}
for pos in LFC_data:
    if CN_d[pos] == 1:
        LFC_n = LFC_data[pos]
        LFC_2n = LFC_data_2n[pos]
        LFC = np.mean([LFC_n,LFC_2n])
        if pos in Pfam_all_d:
            index = Pfam_all_d[pos]
            if index == 0:
                a = 'no_pfam'
            else:
                a = 'yes_pfam'
            if LFC > -1:
                pfam_groups['neutral'][a].append(LFC)
            if -1 > LFC > -3:
                index = Pfam_all_d[pos]
                pfam_groups['weak'][a].append(LFC)
            if -3 >= LFC > -5:
                index = Pfam_all_d[pos]
                pfam_groups['medium'][a].append(LFC)
            if -5 >= LFC:
                index = Pfam_all_d[pos]
                pfam_groups['strong'][a].append(LFC)

outfile = open('pfam_dependent_LFC.txt','w')
outfile.write('group\tattribute\tLFC')
for element in pfam_groups:
    for pfam in pfam_groups[element]:
        for LFC in pfam_groups[element][pfam]:
            outfile.write('\n' + element + '\t' + pfam + '\t' + str(LFC))

