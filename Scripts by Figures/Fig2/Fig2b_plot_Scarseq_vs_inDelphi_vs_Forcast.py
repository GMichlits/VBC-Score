__author__ = 'georg.michlits'

import matplotlib.pyplot as plt
from scipy import stats
import numpy as np

b14k_comp_file = open('b14k_Allen_vs_inDelphi.txt','r')
ScarSeq_d2_inframes = open('inframe_ratio_d2.txt','r')

guide_list = []
colour_list = []
Scar_inframe = []

for line in ScarSeq_d2_inframes:
    col = line.rstrip('\n').split('\t')
    cat = col[0]
    guide = col[1]
    inframe_rat = float(col[2])
    guide_list.append(guide)
    Scar_inframe.append(inframe_rat)
    if cat == 'medium':
        colour_list.append('mediumaquamarine')
    if cat == 'strong':
        colour_list.append('k')
    if cat == 'weak':
        colour_list.append('orange')
    if cat == 'ctrl':
        colour_list.append('silver')

Allen_Shen_dict = {}

i = 0
for line in b14k_comp_file:
    if i > 0:
        col = line.rstrip('\n').split('\t')
        guide = col[0]
        Allen_inframe = float(col[1])
        inDelphi_inframe = float(col[2])
        if guide in guide_list:
            Allen_Shen_dict[guide] = (Allen_inframe, inDelphi_inframe)
    i += 1

Allen_infr_list = []
Shen_infr_list = []
Scar_inframe_plot = []
colour_list_plot = []
for guide in guide_list:
    if guide in Allen_Shen_dict:
        Allen_infr_list.append(Allen_Shen_dict[guide][0])
        Shen_infr_list.append(Allen_Shen_dict[guide][1])
        Scar_inframe_plot.append(Scar_inframe[guide_list.index(guide)])
        colour_list_plot.append(colour_list[guide_list.index(guide)])

print(len(Scar_inframe_plot))
print(len(Allen_infr_list))
print(len(Shen_infr_list))
print(len(colour_list_plot))
#plt.title('Allen')
plt.scatter(Scar_inframe_plot, Allen_infr_list, linewidths=2 , marker='D',edgecolor ='k', s=150, c=colour_list_plot)
gradient, intercept, r_value, p_value, std_err = stats.linregress(Scar_inframe_plot,Allen_infr_list)
print(r_value)
print(p_value)
plt.title('Allen R: ' + str(r_value) + ' p: ' + str(p_value))
mn = 0
mx = 1
x1 = np.linspace(mn,mx,500)
y1 = gradient*x1+intercept
stats.linregress(Scar_inframe_plot, Allen_infr_list)
plt.ylim(0,0.7)
plt.plot(x1, y1, '-r')
plt.savefig('Allen_Forecast_vs_Scarseq.pdf')
plt.close()
plt.close()

#plt.title('Shen')
plt.scatter(Scar_inframe_plot, Shen_infr_list, linewidths=2, marker='D',edgecolor ='k', s=150, c=colour_list_plot)
gradient, intercept, r_value, p_value, std_err = stats.linregress(Scar_inframe_plot,Shen_infr_list)
print(r_value)
print(p_value)
plt.title('Shen R: ' + str(r_value) + ' p: ' + str(p_value))
mn = 0
mx = 1
x1 = np.linspace(mn,mx,500)
y1 = gradient*x1+intercept
stats.linregress(Scar_inframe_plot, Shen_infr_list)
plt.ylim(0,0.7)
plt.plot(x1, y1, '-r')
plt.savefig('Shen_inDelphi_vs_Scarseq.pdf')
plt.close()
plt.close()