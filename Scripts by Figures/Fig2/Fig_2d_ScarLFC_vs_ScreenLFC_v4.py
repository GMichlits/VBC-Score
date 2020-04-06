
__author__ = 'georg.michlits'

dict_file = open('94call_v5_mapv3_indeldict.txt','r')
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import pickle
import math

print('reading guidedepl')
#build guide dictionary
#inframe_ratio_pred[guide][exp][day] = [LFC, p_pos, p_neg]

LFC_screen_dict_n = pickle.load(open('Scarseq_sgRNA_LFCinScreen_n_dict.sav','rb'))
LFC_screen_dict = pickle.load(open('Scarseq_sgRNA_LFCinScreen_dict.sav','rb'))
sgRNA_subgroup_dict = pickle.load(open('Scarseq_sgRNA_subgroups_dict.sav','rb'))

mut_2 = pickle.load(open('mutations_dict_2_top30_relAll_n.sav','rb'))
mut_18 = pickle.load(open('mutations_dict_18_top30_relAll_n.sav','rb'))

##### old stuff not needed anymore but now kept for reference
'''
guide_depl = {}
i = 0
for line in guidedepl_file:
    if i > 0:
        column = line.rstrip('\n').split('\t')
        guide = column[0]
        exp = column[1]
        day = column[2]
        LFC = column[3]
        p_pos = column[4]
        p_neg = column[5]
        if guide not in guide_depl:
            guide_depl[guide] = {}
        if exp not in guide_depl[guide]:
            guide_depl[guide][exp] = {}
        if day not in guide_depl[guide][exp]:
            guide_depl[guide][exp][day] = [LFC, p_pos, p_neg]
    i += 1
### determine scarSeq LFC
mut_2 = pickle.load(open('mutations_dict_2_top30_relAll_n.sav','rb'))
mut_18 = pickle.load(open('mutations_dict_18_top30_relAll_n.sav','rb'))

exclude_list = ['Line1_1','Line1_2','Doxl2_6','ROSA26_1','Rhox2f_3','Rhox2f_1','ROSA26_2','ROSA26_3','Cdk1_2','Olfr10_1','Olfr52_1','Olfr44_1','Traip_4','Traip_7']
weak = ['Srf_2','Rhox2f_1','Traip_2','Eif3b_5','Hscb_3','Dcaf13_2','Nup214_4','Traip_3','Clns1a_1','Jak1_2','Traip_4','Taf8_5','Sub1_2','Chek1_1']
medium = ['Jak1_5','Snapc5_3','Eif3a_5','Taf8_4','Cdk1_2','Sub1_1','Dcaf13_4','Clns1a_5','Rpl35a_4','Rpl35a_6','Snapc5_1','Nup214_2','Eif3b_4']
Strong = ['Ruvbl2_1','Rpl13a_1','Psmd2_4','Nelfb_6','Rpl13a_4','Kif18b_2','Rad51_4','Psma1_4','Kif11_5','Kif18b_4','Ggps1_3','Eif3a_4','Rad51_1','Cdk1_4','Psma2_2','Nelfb_2','Spc25_6','Socs3_2','Traip_1','Rrm1_4','Ggps1_2','Socs3_4','Srf_3','Kif11_1','Ruvbl2_2','Ncbp1_5','Psma2_3','Txn1_2','Rrm1_2','Txn1_1','Stat3_4','Traip_6','Stat3_2','Rack1_1','Spc25_2','Psma1_3','Ncbp1_2','Rack1_2','Chek1_5','Psmd2_3','Hscb_1','Rbx1_1']
CTRL = ['ROSA26_1','ROSA26_2','ROSA26_3','Olfr10_1','Olfr52_1','Olfr44_1']
exclude_list = ['Line1_1','Line1_2','ROSA26_1','Rhox2f_3','Rhox2f_1','ROSA26_2','ROSA26_3','Olfr10_1','Olfr52_1','Olfr44_1','Traip_7']
'''
##### old stuff not needed anymore but now kept for reference

LFC_scar = {}
for guide in mut_18:
    LFC_scar[guide] = 0
    day18_del_frac = 0
    day18_ins_frac = 0
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
                day18_del_frac += np.mean(np.array([repl_1_value,repl_2_value,repl_3_value])) #3 replicates
        if type == 'ins':
            if -21 < int(pos_start) < 21:
                x = 1*int(dist)
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
                day18_ins_frac += np.mean(np.array([repl_1_value,repl_2_value,repl_3_value])) #3 replicates
    print(guide + '\t' + str(day18_del_frac+day18_ins_frac))
    frac_left = day18_del_frac+day18_ins_frac
    if frac_left == 0:
        frac_left = 0.01
    LFC_scar[guide] = math.log2(frac_left)

x = []
y = []
c_list = []
# guides with name CTRL_.. are anyway not in predicted_inframeratio guide dict so they are not part of the graph.
for guide in LFC_scar:
    if not LFC_scar[guide] == []:
        if guide in LFC_screen_dict:
            inframe_bit_value = np.median(LFC_scar[guide])
            depl = LFC_screen_dict[guide]
            y.append(float(depl))
            x.append(float(inframe_bit_value))
            subgroup = sgRNA_subgroup_dict[guide]
            if subgroup == 'weak':
                c = 'orange'
            elif subgroup == 'medium':
                c = 'mediumaquamarine'
            elif subgroup == 'strong':
                c = 'black'
            c_list.append(c)

name = 'LFC_Screen_vs_LFC_ScarSeq'
plt.title(name)
plt.figure(figsize=(7, 7))
plt.scatter(x, y, c=c_list)
plt.ylabel('LFC_screen')
plt.xlabel('LFC_scarSeq')
a = np.array(x)
b = np.array(y)
#gradient, intercept, r_value, p_value, std_err = stats.linregress(a,b)
#plt.text(-10, .95, 'slope '+ str(gradient))
#plt.text(-10, .9, 'Rsqr '+ str(r_value))
#mn=np.min(a)
#mx=np.max(b)
#x1=np.linspace(mn,mx,500)
#y1=gradient*x1+intercept
#plt.plot(a,b,'ob')
#plt.plot(x1,y1,'-r')

#plt.show()
plt.savefig('v5.3_data/' + name + 'v4.pdf')
