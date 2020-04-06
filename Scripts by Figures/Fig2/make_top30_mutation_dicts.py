
__author__ = 'georg.michlits'

dict_file = open('94call_v5_mapv3_indeldict.txt','r')
import numpy as np
import pickle


#build guide dictionary
#inframe_ratio_pred[guide][exp][day] = [LFC, p_pos, p_neg]

exp_type_plt = 'n'
top_x = 30
pos_range = 21

print('assemble indeldict')
#indelDICT[guide][repl][type][d][mutation]=count
INDEL_DICT ={}
i=0
dict_file.seek(0)
for line in dict_file:
    if i > 0:
        column = line.rstrip('\n').split('\t')
        guide = column[0]
        repl = column[1]
        exp_type = column[2]
        d = column[3].rstrip(' ')
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

mutation_dict_18 = {}
mutation_dict_2 = {}
total_count = {}
wt_count = {}
total_count_18 = {}
wt_count_18 = {}
for guide in INDEL_DICT:
#if True:
    #guide = 'Dcaf13_4'
    x_list = []
    x_label = []
    y_list = []
    yerr_list = []
    c_list = []
    Freq_dict = {}
    mutation_dict_2[guide] = {}
    for repl in INDEL_DICT[guide]:
        #for d in INDEL_DICT[guide][repl][exp][day]:
        if 'day2' in INDEL_DICT[guide][repl][exp_type_plt]:
            #mutation_dict_2[guide][repl][mutation] = frequ
            x=0
            total_count[guide]={}
            wt_count[guide] = {}
            total_count[guide][repl]=0
            wt_count[guide][repl]=0
            for mutation_item in sorted(INDEL_DICT[guide][repl][exp_type_plt]['day2'].items(), key=lambda t: float(t[1]), reverse=True):
                mutation = mutation_item[0]
                count = float(mutation_item[1])
                if mutation.split(':')[0] == 'ins' or mutation.split(':')[0] == 'del':
                    # select only guides with pos +- 21 bp of cut size this exclude e.g. del in primers that would giv del:1:-132 mutations and similar
                    pos = int(mutation.split(':')[2])
                    if -1*pos_range <= pos <= pos_range:
                        if mutation not in mutation_dict_2[guide]:
                            mutation_dict_2[guide][mutation] = {}
                        total_count[guide][repl] += count
                        mutation_dict_2[guide][mutation][repl] = count
                if mutation.split(':')[0] == 'wt':
                    wt_count[guide][repl] += count

            x = 0
            x_list = []
            x_label = []
            y_list = []
            c_list = []
            for mutation in mutation_dict_2[guide]:
                if repl in mutation_dict_2[guide][mutation]:
                    count = float(mutation_dict_2[guide][mutation][repl])
                else:
                    count = 0
                    mutation_dict_2[guide][mutation][repl] = 0
                #now change count to frequencies
                total_count_v = (wt_count[guide][repl]+total_count[guide][repl])
                if total_count_v == 0:
                    mutation_dict_2[guide][mutation][repl] = 0
                else:
                    mutation_dict_2[guide][mutation][repl] = count/total_count_v
    for mutation in mutation_dict_2[guide]:
        values_mean = []
        for repl in mutation_dict_2[guide][mutation]:
            values_mean.append(mutation_dict_2[guide][mutation][repl])
        freq_mean = np.mean(values_mean)
        freq_std = np.std(values_mean)
        Freq_dict[mutation] = (freq_mean,freq_std)
    x = 0
    len_mutations = len(Freq_dict)
    for mut_item in sorted(Freq_dict.items(), key=lambda t: t[1][0], reverse=True):
        mutation = mut_item[0]
        y = mut_item[1][0]
        yerr = mut_item[1][1]
        if x < min(top_x, len_mutations):
            x += 1
            x_list.append(x)
            y_list.append(y)
            yerr_list.append(yerr)
            if int(mutation.split(':')[1])%3 == 0:
                c = 'g'
            else:
                c = 'r'
            c_list.append(c)
            mutation = mutation.replace('del','d')
            mutation = mutation.replace('ins','i')
            x_label.append(mutation)

    #repeat generating negative data for day 18 (use the most freuqent indels found at day2)
    y_list_18 = []
    yerr_list_18 = []
    Freq_dict_18 = {}
    mutation_dict_18[guide] = {}
    for repl in INDEL_DICT[guide]:
        #for d in INDEL_DICT[guide][repl][exp][day]:
        if 'day2' in INDEL_DICT[guide][repl][exp_type_plt]:
            #mutation_dict_2[guide][repl][mutation] = frequ
            x = 0
            total_count_18[guide] = {}
            wt_count_18[guide] = {}
            total_count_18[guide][repl] = 0
            wt_count_18[guide][repl] = 0
            for mutation_item in sorted(INDEL_DICT[guide][repl][exp_type_plt]['day2'].items(), key=lambda t: float(t[1]), reverse=True):
                mutation = mutation_item[0]
                if mutation in INDEL_DICT[guide][repl][exp_type_plt]['day18']:
                    count = float(INDEL_DICT[guide][repl][exp_type_plt]['day18'][mutation])
                else:
                    count = 0
                if mutation.split(':')[0] == 'ins' or mutation.split(':')[0] == 'del':
                    # select only guides with pos +- 21 bp of cut size this exclude e.g. del in primers that would giv del:1:-132 mutations and similar
                    pos = int(mutation.split(':')[2])
                    if -1*pos_range <= pos <= pos_range:
                        if mutation not in mutation_dict_18[guide]:
                            mutation_dict_18[guide][mutation] = {}
                        total_count_18[guide][repl] += count
                        mutation_dict_18[guide][mutation][repl] = count
                if mutation.split(':')[0] == 'wt':
                    wt_count_18[guide][repl] += count

            for mutation in mutation_dict_18[guide]:
                if repl in mutation_dict_18[guide][mutation]:
                    count = float(mutation_dict_18[guide][mutation][repl])
                else:
                    count = 0
                    mutation_dict_18[guide][mutation][repl] = 0
                #now change count to frequencies    # the factor 3.6 is emperically determined to normalize for ROSA, and Olfactory receptor enrichment of silenced guides - no DNA strand breaks
                # used either wt_count or total_count as normalization - total count is more robust but makes less efficient guides less important
                total_count_v = (wt_count_18[guide][repl]+total_count_18[guide][repl])/3.6
                if total_count_v == 0:
                    mutation_dict_18[guide][mutation][repl] = 0
                else:
                    mutation_dict_18[guide][mutation][repl] = count/total_count_v
    for mutation in mutation_dict_18[guide]:
        values_mean = []
        for repl in mutation_dict_18[guide][mutation]:
            values_mean.append(mutation_dict_18[guide][mutation][repl])
        freq_mean = np.mean(values_mean)
        freq_std = np.std(values_mean)
        Freq_dict_18[mutation] = (freq_mean,freq_std)

    len_mutations = len(Freq_dict)
    x = 0
    for mut_item in sorted(Freq_dict.items(), key=lambda t: t[1][0], reverse=True):
        mutation = mut_item[0]
        data_18 = Freq_dict_18[mutation]
        y_18 = data_18[0] * -1
        yerr_18 = data_18[1]
        if x < min(top_x, len_mutations):
            x += 1
            y_list_18.append(y_18)
            yerr_list_18.append(yerr_18)

pickle.dump(mutation_dict_2,open('mutations_dict_2_top30_relwt_n.sav','wb'))
pickle.dump(mutation_dict_18,open('mutations_dict_18_top30_relwt_n.sav','wb'))


