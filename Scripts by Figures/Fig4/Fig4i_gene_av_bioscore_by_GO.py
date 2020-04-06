__author__ = 'georg.michlits'

import numpy as np
import statistics as st
import pickle
Gene_to_20nt_to_pos = pickle.load(open('../../Gen_properties/Ms/property_sav/Gene_to_20nt_to_pos.sav', 'rb'))
AA_cons_d = pickle.load(open('../../Gen_properties/Ms/property_sav/Jun_Bioscore_3Nov_d.sav', 'rb'))


Ms_cons = open('Msgene_Bioscore.txt','r')
GO_descr = open('goTermsDescription.txt','r')
gene_Go = open('gene_association.mgi.txt','r')
CEG2_gene_set = open('CEG2_mouse.txt','r')
outfile = open('Bioscore_vs_mean_std.txt','w')

i = 0
Essentials = {}
for line in CEG2_gene_set:
    gene = line.rstrip('\n')
    Essentials[gene] = 0
    i += 1

Ms_con_dic = {}
i = 0
Ms_biosc_max_gd = {}
for line in Ms_cons:
    if i > 0:
        column = line.rstrip('\n').split('\t')
        gene = column[0]
        cons = column[1]
        Ms_biosc_max_gd[gene] = float(column[2])
        if not cons == 'nan':
            Ms_con_dic[gene] = float(cons)
    i += 1

GO_gene_dict = {}
for line in gene_Go:
    if not line.startswith('!'):
        column = line.rstrip('\n').split('\t')
        gene = column[2]
        GO_ID = column[4]
        if GO_ID not in GO_gene_dict:
            GO_gene_dict[GO_ID] = []
        if gene not in GO_gene_dict[GO_ID]:
            GO_gene_dict[GO_ID].append(gene)

outfile.write('GO_ID\tGO_term\tGO_class\tn_genes\tmean_bioscore\tmaxbioscore\twithin_gene_Bioscore_std\tstd_cons\tsem_cons\tgenelist_max100')
count_goterms = 0

for line in GO_descr:
    count_goterms += 1
    print(count_goterms)
    column = line.rstrip('\n').split(':')
    GO_ID = column[0] + ':' + column[1]
    if '*' in column[2]:
        GO_term = column[2].split('*')[0]
        GO_class = column[2].split('*')[1]
    else:
        GO_term = column[2]
        GO_class = 'NA'
    if GO_ID in GO_gene_dict:
        gene_list = GO_gene_dict[GO_ID]
        number_of_genes = len(gene_list)
        cons_list = []
        STD_AAcons_list = []
        pass_genelist = []
        max_Bioscore_list = []
        k = 0
        for gene in gene_list:
            if gene in Ms_con_dic:
                max_Bioscore_list.append(Ms_biosc_max_gd[gene])
                AA_list = []
                if gene in Gene_to_20nt_to_pos:
                    for _20nt in Gene_to_20nt_to_pos[gene]:
                        pos = Gene_to_20nt_to_pos[gene][_20nt]
                        if pos in AA_cons_d:
                            AA_cons = AA_cons_d[pos]
                            if AA_cons > 0:
                                AA_list.append(AA_cons)
                if len(pass_genelist) < 100:
                    pass_genelist.append(gene)
                if len(AA_list) > 10:
                    STD_AAcons = st.stdev(AA_list)
                    STD_AAcons_list.append(STD_AAcons)
                k+=1
                cons_list.append(float(Ms_con_dic[gene]))
        if k >= 15:
            mean = st.mean(cons_list)
            within_gene_std_AAcons = st.mean(STD_AAcons_list)
            max_Bioscore_GO = st.mean(max_Bioscore_list)
            std = st.stdev(cons_list)
            sem = st.stdev(cons_list)/np.sqrt(len(cons_list))
            outfile.write('\n' + str(GO_ID))
            outfile.write('\t' + str(GO_term))
            outfile.write('\t' + str(GO_class))
            outfile.write('\t' + str(k))
            outfile.write('\t' + str(mean))
            outfile.write('\t' + str(max_Bioscore_GO))
            outfile.write('\t' + str(within_gene_std_AAcons))
            outfile.write('\t' + str(std))
            outfile.write('\t' + str(sem))
            outfile.write('\t' + str(pass_genelist))


cons_list = []
pass_genelist = []
STD_AAcons_list = []
max_Bioscore_list = []
k = 0
for gene in Essentials:
    if gene in Ms_con_dic:
        max_Bioscore_list.append(Ms_biosc_max_gd[gene])
        AA_list = []
        hiAA = 0
        if gene in Gene_to_20nt_to_pos:
            for _20nt in Gene_to_20nt_to_pos[gene]:
                pos = Gene_to_20nt_to_pos[gene][_20nt]
                if pos in AA_cons_d:
                    AA_cons = AA_cons_d[pos]
                    if AA_cons > 0:
                        AA_list.append(AA_cons)
        if len(pass_genelist) < 100:
            pass_genelist.append(gene)
        if len(AA_list) > 10:
            STD_AAcons = st.stdev(AA_list)
            for element in AA_list:
                if element > 0.5:
                    hiAA += 1
            STD_AAcons_list.append(STD_AAcons)
        k+=1
        cons_list.append(float(Ms_con_dic[gene]))
if k >= 5:
    mean = st.mean(cons_list)
    within_gene_std_AAcons = st.mean(STD_AAcons_list)
    std = st.stdev(cons_list)
    max_Bioscore_GO = st.mean(max_Bioscore_list)
    sem = st.stdev(cons_list)/np.sqrt(len(cons_list))
    outfile.write('\n' + str('CEG2_Hart'))
    outfile.write('\t' + str('CEG2_Hart'))
    outfile.write('\t' + str('CEG2_Hart'))
    outfile.write('\t' + str(k))
    outfile.write('\t' + str(mean))
    outfile.write('\t' + str(max_Bioscore_GO))
    outfile.write('\t' + str(within_gene_std_AAcons))
    outfile.write('\t' + str(std))
    outfile.write('\t' + str(sem))
    outfile.write('\t' + str(pass_genelist))


