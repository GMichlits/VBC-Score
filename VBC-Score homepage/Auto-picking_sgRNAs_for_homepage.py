__author__ = 'georg.michlits'

import pickle
import re
import numpy as np

species_list = ['Ms','Hm']#, 'Hm']
#species_list = ['mm10', 'ce11','rn6', 'xenTro9', 'dm6']
#species_list = ['xenTro9']
#species_list = ['ce11']
#species_list = ['rn6', 'xenTro9', 'dm6']

# SETTINGS:
topX_sgRNA_autopickmode = 50

#### PENALTIES for position of sgRNAs can be used to avoid targeting the very same region in the protein multiple times
# it pick the best sgRNA first for the next picks it gives penalites to nearby sgRNAs (within a window of x % CDS)
# region of X% where sgRNAs get penalties if better guide targeting the same region was already picked
region_window = 0.05
penality_2ndsgRNA = 1 #with 1 the feature is deactivated
penality_3rdsgRNA = 1 #with 1 the feature is deactivated
# penality for sgRNAs if only 1bp away. E.g. To avoid to target a region with -GGGG- Sequence with 3 separate sgRNAs all
# 1 bp shifted to one another
#proximity_lessthanonebp_penality = 1 #with 1 the feature#  is deactivates
proximity_lessthantwobp_penality = 0.5

# which length of sgRNA is preferred? 20nt natural G is most active shorter sgRNAS might avoid off-targets at the cost of reducing activity
# G-overhang also dampens sgRNA activity
penalty_19nt = 1.05
penalty_18nt = 1.0
penalty_20nt = 1.1
penalty_Gplus20nt = 0.95

# poly T streches act as PolIII terminator and should be avoided
penalty_5Tstretch = 0.01

#avoiding sgRNAS with OFF-targets (only looking at annotated Genes - coding DNA Sequence)
penality_OT95over0 = 0.3
penalty_OT90over10 = 0.05
penalty_OT90_1to3 = 0.6
penalty_OT90_3to6 = 0.3
penalty_OT90_6to10 = 0.1

#for most genes multiple isoforms exists - this gives a penalty to sgRNAs for targeting Isoform-specific exons
penality_noncommon_CDS_isoforms = 0.1

# a penality for sgRNAs that have a perfect match (entire genome including introns and intergene space)
# (sgRNAs likely to deplete because of copy number effect)
penalty_2CN = 0.9
penalty_over3CN = 0.75

def pos_to_cutsite(pos):
    orientation = pos.split('(')[1][0]
    if orientation == '+':
        start_pos = int(pos.split(':')[1].split('-')[0])
        cutsite = start_pos + 21
    elif orientation == '-':
        start_pos = int(pos.split('-')[1].split('(')[0])
        cutsite = start_pos - 21
    else:
        print('ERROR_guide has no orientation' )
    return cutsite
def revcomp(dna, reverse=True, complement=True):
    DNA = dna.upper()
    a = 'ACGTRYKMSWBDHVN'
    b = 'TGCAYRMKSWVHDBN'
    dict = {a[i]:b[i] for i in range (15)}
    if reverse:
        DNA = reversed(DNA)
    if complement:
        result = [dict[i] for i in DNA]
    else: result = [i for i in DNA]
    return ''.join(result)

for Species in species_list:
    print(Species)
    basic_property_filename = Species + '_mergeData_all_ext80bp.txt'
    Basic_property_file = open('../../Gen_properties/' + Species + '/' + basic_property_filename, 'r')

    outfile_top6 = open('top_6/' + Species + '_' + str(topX_sgRNA_autopickmode) + '_sgRNAs_Dec12.txt','w')

    Gene_d = pickle.load(open('../../Gen_properties/' + Species + '/property_sav/Gene_d.sav', 'rb'))
    distance_d = pickle.load(open('../../Gen_properties/' + Species + '/property_sav/distance_d.sav', 'rb'))
    CN_d = pickle.load(open('../../Gen_properties/' + Species + '/property_sav/CN_d.sav', 'rb'))
    combo_sgRNA_d = pickle.load(open('../../Gen_properties/' + Species + '/property_sav/combo_sgRNA_d.sav', 'rb'))
    inframe_InDelphi = pickle.load(open('../../Gen_properties/' + Species + '/property_sav/inDelphi_infr_d.sav', 'rb'))
    VBC_d = pickle.load(open('../../Gen_properties/' + Species + '/property_sav/Jul14_VBC_d.sav', 'rb'))
    Bioscore_d = pickle.load(open('../../Gen_properties/' + Species + '/property_sav/Jul14_Bioscore_d.sav', 'rb'))
    alias_gene_names_gd = pickle.load(open('../../Gen_data/VBC-score compatibility/' + Species + 'alias_gd.sav', 'rb'))

    print('extract AA-cut and nt30 from properties File')
    nt30_d = {}
    pep_d = {}
    res_d = {}
    common_CDS_d = {}
    count_guides = 0
    multi_Ts_d = {}
    for i,line in enumerate(Basic_property_file):
        if i > 0:
            col = line.rstrip('\n').split('\t')
            pos = col[1]
            pep = col[26]
            common_CDS = col[16]
            common_CDS_d[pos] = common_CDS
            pep_d[pos] = pep
            nt30 = col[3]
            nt30_d[pos] = nt30.upper()
            if not re.search(r"T{4,25}", nt30[4:27]) == None:
                multi_Ts_d[pos] = 1
            else:
                multi_Ts_d[pos] = 0
            if 'GAAGAC' in nt30[4:27] or 'CTTCTG' in nt30[4:27] or 'CGTCTC' in nt30[4:27] or 'GCAGAG' in nt30[4:27]:
                res_d[pos] = 'warning BsbI or BsmBI restriction site in sgRNA seq'
            elif multi_Ts_d[pos] == 1:
                res_d[pos] = 'multi-T-strech in sgRNA - 4 or more PolIII terminator'
            else:
                res_d[pos] = '-'
    Basic_property_file.close()

    print('read in Off-target prediction')
    if Species == 'Hm':
        off_tartet_species = 'hg38'
    elif Species == 'Ms':
        off_tartet_species = 'mm10'
    else:
        off_tartet_species == Species
    OFFtarget_filename = off_tartet_species + '_offtarget.txt'
    OFFtarget_file = open('../../VBC_cross_species/' + OFFtarget_filename)
    OT90_d = {}
    OT95_d = {}
    OT90_offlist_d = {}
    OT95_offlist_d = {}
    for line in OFFtarget_file:
        col = line.rstrip('\n').split('\t')
        pos = col[0]
        OT95 = int(col[1])
        OT90 = int(col[2])
        OT95_pos_list = col[5].split('|')
        OT90_pos_list = col[6].split('|')
        OT95_pos_list.pop(OT95_pos_list.index('')) #Thomas adds a | also to the last enrty therefore i pop the last empyt enrty
        OT90_pos_list.pop(OT90_pos_list.index(''))
        OT90_d[pos] = OT90
        OT95_d[pos] = OT95
        OT90_offlist_d[pos] = OT90_pos_list
        OT95_offlist_d[pos] = OT95_pos_list
    OFFtarget_file.close()

    print('get AA-targets from new properties')
    if Species == 'Hm':
        off_tartet_species = 'hg38'
    elif Species == 'Ms':
        off_tartet_species = 'mm10'
    new_property_file = open('../../VBC_cross_species/' + off_tartet_species + '_cut.propertie.new2.txt', 'r')
    for i,line in enumerate(new_property_file):
        if i > 0:
            col = line.rstrip('\n').split('\t')
            pos = col[1]
            pep = col[6]
            pep_d[pos] = pep
    new_property_file.close()

    #assign by cross-checking GeneRef gtf file
    #exon_dict[chromosome ID][gene_name or gene_ID]=[start,stop,exon]
        # if guide start or guide stop between annotated exon assign exon = else write NA.
    print('read in exon assignements')
    if Species == 'Hm':
        gtf_file = open('../../Gen_properties/Hm/hg38_refGene.gtf','r')
    elif Species == 'Ms':
        gtf_file = open('../../Gen_properties/Ms/mm10_refGene.gtf','r')
    exon_dict = {}
    for line in gtf_file:
        col = line.rstrip('/n').split('\t')
        chr = col[0]
        if col[2] == 'exon':
            start = int(col[3])
            stop = int(col[4])
            exon = col[8].split(';')[2].replace('_number' ,'').replace('"', '')
            NM_ID = col[8].split(';')[0].replace('gene_id ' ,'').replace('"', '')
            gene_ID = col[8].split(';')[0].replace('gene_id ' ,'').replace('"', '')
            gene_name = col[8].split(';')[4].replace('gene_name ' ,'').replace('"', '')
            if not chr in exon_dict:
                exon_dict[chr] = {}
            if not gene_name in exon_dict[chr]:
                exon_dict[chr][gene_name] = []
            if not gene_ID in exon_dict[chr]:
                exon_dict[chr][gene_ID] = []
            if gene_ID == gene_name:
                exon_dict[chr][gene_ID].append([start,stop,exon])
            else:
                exon_dict[chr][gene_ID].append([start,stop,exon])
                exon_dict[chr][gene_name].append([start,stop,exon])
            # i do this twice for those cases were gene id and genename are equal chr and exon position are filters anyway
    gtf_file.close()

    print('adjust VBC-score for top6 selection "VBC-score-picking"')
    sgRNAs = {}
    natural_g_dict = {18:{},19:{},20:{}}
    for pos in nt30_d:
        chromosome_name = pos.split(':')[0]
        if not 'Un' in pos.split(':')[0] and not 'N' in nt30_d[pos] and not '' == nt30_d[pos]:
            OT90 = OT90_d[pos]
            OT95 = OT95_d[pos]
            CN = CN_d[pos]
            #print(pos + '\t' + nt30_d[pos] +'\t' + Gene_d[pos])
            VBC = VBC_d[pos]
            Ts = multi_Ts_d[pos]
            common_CDS = common_CDS_d[pos]
            seq = nt30_d[pos]
            # determine if it is possible to make a guide that starts with natural G (for length 20bp 19 bp or 18bp):
            if seq[4] == 'G':
                naturalG_20 = 1
                natural_g_dict[20][pos] = 0
            else:
                naturalG_20 = 0
            if seq[5] == 'G':
                naturalG_19 = 1
                natural_g_dict[19][pos] = 0
            else:
                naturalG_19 = 0
            if seq[6] == 'G':
                naturalG_18 = 1
                natural_g_dict[18][pos] = 0
            else:
                naturalG_18 = 0
            VBC_score_picking = VBC
            # penality for copy number (number of perfect target sites in the genome):
            if CN == 2:
                VBC_score_picking = VBC_score_picking * penalty_2CN
            elif CN >= 3:
                VBC_score_picking = VBC_score_picking * penalty_over3CN
            # penality for Off-target (OT95 is an estimate of very likely off-target sites in the genome):
            if OT90 > 0:
                VBC_score_picking = VBC_score_picking * penality_OT95over0
            # penality for Off-target (OT90 is an estimate of likely off-target sites in the genome):
            if 0 < OT90 <=3:
                VBC_score_picking = VBC_score_picking * penalty_OT90_1to3
            elif 3 < OT90 <= 6:
                VBC_score_picking = VBC_score_picking * penalty_OT90_3to6
            elif 6 < OT90 <= 10:
                VBC_score_picking = VBC_score_picking * penalty_OT90_6to10
            elif OT90 > 10:
                VBC_score_picking = VBC_score_picking * penalty_OT90over10
            # penality for TTTTT strech (5T in a row act as PolIII terminator sequence and should be avoided):
            if Ts == 1:
                VBC_score_picking = VBC_score_picking * penalty_5Tstretch
            ######  This is only required if non-common exons are considered as well
            if common_CDS == 'FALSE':
                VBC_score_picking = VBC_score_picking * penality_noncommon_CDS_isoforms
            if naturalG_19 == 1:
                VBC_score_picking = VBC_score_picking * penalty_19nt
            elif naturalG_18 == 1:
                VBC_score_picking = VBC_score_picking * penalty_18nt
            elif naturalG_20 == 1:
                VBC_score_picking = VBC_score_picking * penalty_20nt
            else:
                VBC_score_picking = VBC_score_picking * penalty_Gplus20nt
            genename = Gene_d[pos]
            if not genename in sgRNAs:
                sgRNAs[genename] = {}
            sgRNAs[genename][pos] = VBC_score_picking
            #print(genename + '\t' + str(VBC) + '\t' +str(VBC_score_picking) + '\t' + str(OT90) + '\t' + str(CN) + '\t' + str(Ts))

    outfile_top6.write('gene' + '\t' + 'sgRNA' + '\t' + 'Position' + '\t' + 'Exon number refGene' + '\t' + 'Auto-pick top sgRNAs' + '\t' + 'VBC score' + '\t'
                       + 'Distance TSS=0 stop=1' + '\t' +'sgRNA activity 0=bad 1=good' + '\t' + 'Frameshift ratio inDelphi' + '\t' +
                        'Bioscore 0=bad 1=good' + '\t' + 'Likely off-targets OT>90' + '\t' +
                        'Cutting site amino acids' + '\t' +
                        'Very likely off-targets OT>95' + '\t' +
                        'Very likely off-target genes' + '\t' + 'Likely off-targets OT>90' + '\t' +
                        'Likely off-target genes' + '\t' + 'Perfect match sites in genome' + '\t' +
                        'Recommended sgRNA-length' + '\t'
                        'GeCKO cloning FW_oligo' + '\t' + 'GeCKO cloning RV_oligo' + '\t' + 'Additional properties' + '\t' +
                        'alias gene names')

    #outfile_top6.write('genename' + '\t' + 'position' + '\t' + 'automated selection "top 6"'+ '\t' + '20nt+NGG' + '\t' + 'exon number refGene' + '\t'
    #                    'recommended sgRNA-length' + '\t' + 'cutting site amino acids' + '\t' +
    #                    'distance TSS=0 to stop=1' + '\t' + 'VBC score 0=bad_1=good' + '\t' +
    #                    'sgRNA_activity 0=bad 1=good' + '\t' + 'frameshift prediction inDelphi' + '\t' +
    #                    'Bioscore 0=bad 1=good' + '\t' + 'very likely off-targets OT>95' + '\t' +
    #                    'very likely off-target genes' + '\t' + 'likely off-targets OT>90' + '\t' +
    #                    'likely off-target genes' + '\t' + 'perfect match sites in genome' + '\t' +
    #                    'GeCKO cloning FW_oligo' + '\t' + 'GeCKO cloning RV_oligo' + '\t' + 'CAUTION BbSI or BsmBI site in sgRNA')

    print('run selection process pick ' + str(topX_sgRNA_autopickmode) + ' sgRNA per gene')
    # here we run the selection for every gene 6 times (to get 6 guides) one by one we modify the vbc_score_picking values
    # depending on proximity to regions that have been target with previously every selected guide varies the vbc_score_picking_values.
    top6_d = {}
    for gene in sorted(sgRNAs):
        #hitzones are entered as lists [start,stop]
        #print(gene)
        hitzone = []
        for i in range(topX_sgRNA_autopickmode):
            for item in sorted(sgRNAs[gene].items(), key = lambda t: t[1], reverse=True):
                pos = item[0]
                guide_selection_rank = i+1
                top6_d[pos] = guide_selection_rank
                cutsite_picked_sgRNA = pos_to_cutsite(pos)
                orientation_picked_sgRNA = pos.split('(')[1][0]
                VBC_score_picking = item[1]
                dist = distance_d[pos]
                hit_type = '1st'
                # check if it is a 'second hit within a 10% range of a previous hit zone in the protein
                for element in hitzone:
                    if element[0] < dist < element[1]:
                        hit_type = '2nd'
                        # extend the size of the zone depending on the new hit
                        new_zone = [min(element[0],dist-region_window/2),max(element[1],dist+region_window/2)]
                        hitzone[hitzone.index(element)] = new_zone
                        break
                # remove selected sgRNA from dictionary
                sgRNAs[gene].pop(pos)
                ####################################################
                ##########     get guide properties and write result
                ####################################################
                #print(pos + '\t' + str(dist))
                if pos in natural_g_dict[19]:
                    guide_sequence = nt30_d[pos][5:24]
                    oligotype = '19nt'
                elif pos in natural_g_dict[20]:
                    guide_sequence = nt30_d[pos][4:24]
                    oligotype = '20nt'
                elif pos in natural_g_dict[18]:
                    guide_sequence = nt30_d[pos][6:24]
                    oligotype = '18nt'
                else:
                    guide_sequence = 'G' + nt30_d[pos][4:24]
                    oligotype = 'G + 20nt'

                FW_oligo = 'CACC' + guide_sequence
                RV_oligo = revcomp(guide_sequence + 'GTTT')


                # from OT dictionaries retrieve information nt seq and genename of OT95 (very likely) or OT90 likely off targets
                OT95_offtargets = '0'
                for n_element,element_pos in enumerate(OT95_offlist_d[pos]):
                    if element_pos in Gene_d:
                        genename = Gene_d[element_pos]
                    else:
                        genename = 'N.A.'
                    if element_pos in nt30_d:
                        nt_23 = nt30_d[element_pos][4:27]
                    else:
                        nt_23 = element_pos
                    if n_element == 0:
                        OT95_offtargets = genename + '|' + nt_23
                    else:
                        OT95_offtargets += ',' + genename + '|' + nt_23
                OT90_offtargets = '0'
                for n_element,element_pos in enumerate(OT90_offlist_d[pos]):
                    if element_pos in Gene_d:
                        genename = Gene_d[element_pos]
                    else:
                        genename = 'N.A.'
                    if element_pos in nt30_d:
                        nt_23 = nt30_d[element_pos][4:27]
                    else:
                        nt_23 = element_pos
                    if n_element == 0:
                        OT90_offtargets = genename + '|' + nt_23
                    else:
                        OT90_offtargets += ',' + genename + '|' + nt_23

                # get exon number fro exon dict
                # exon_dict[chromosome ID][gene_name or gene_ID]=[start,stop,exon]
                chr = pos.split(':')[0]
                start = int(pos.split(':')[1].split('-')[0])
                stop = int(pos.split(':')[1].split('-')[1].split('(')[0])
                exon = 'NA'
                for test_set in exon_dict[chr][gene]:
                    test_start = test_set[0]
                    test_stop = test_set[1]
                    test_exon = test_set[2]
                    if test_start < start < test_stop or test_start < stop < test_stop:
                        exon = test_exon
                        break #breaks the loop for test in exon... because exon was already found
                    else:
                        exon = 'isoform specific exon'
                # exon = NA or if assigned 'exon 3' (or similar)
                if not pos in OT95_d:
                    OT95_d[pos] = 'NA'
                if not pos in OT90_d:
                    OT90_d[pos] = 'NA'
                if not pos in CN_d:
                    CN_d[pos] = 'NA'
                if not pos in pep_d:
                    pep_d[pos] = 'NA'
                #print(gene + '\t' + nt30_d[pos][4:27] + '\t' + str(OT90_offtargets))
                if gene in alias_gene_names_gd:
                    alias = str(alias_gene_names_gd[gene])
                else:
                    alias = ''
                if not pos in inframe_InDelphi:
                    outfile_top6.write('\n' + gene + '\t' + nt30_d[pos][4:27] + '\t'
                                       + pos + '\t'
                                       + exon + '\t'
                                       + str(guide_selection_rank) + '\t'
                                       + str(round(VBC_d[pos],3)) + '\t'
                                       + str(round(distance_d[pos],3)) + '\t'
                                       + str(round(combo_sgRNA_d[pos],3)) + '\t'
                                       + str('NA') + '\t'
                                       + str(round(Bioscore_d[pos],3)) + '\t'
                                       + str(OT90_d[pos]) + '\t'
                                       + pep_d[pos] + '\t'
                                       + str(OT95_d[pos]) + '\t' + OT95_offtargets + '\t'
                                       + str(OT90_d[pos]-OT95_d[pos]) + '\t' + OT90_offtargets + '\t'
                                       + str(CN_d[pos]) + '\t'
                                       + oligotype + '\t' + FW_oligo + '\t' + RV_oligo + '\t' + res_d[pos] + '\t'
                                       + alias)
                else:
                    outfile_top6.write('\n' + gene + '\t' + nt30_d[pos][4:27] + '\t'
                                       + pos + '\t'
                                       + exon + '\t'
                                       + str(guide_selection_rank) + '\t'
                                       + str(round(VBC_d[pos],3)) + '\t'
                                       + str(round(distance_d[pos],3)) + '\t'
                                       + str(round(combo_sgRNA_d[pos],3)) + '\t'
                                       + str(round(1-inframe_InDelphi[pos],3)) + '\t'
                                       + str(round(Bioscore_d[pos],3)) + '\t'
                                       + str(OT90_d[pos]) + '\t'
                                       + pep_d[pos] + '\t'
                                       + str(OT95_d[pos]) + '\t' + OT95_offtargets + '\t'
                                       + str(OT90_d[pos]-OT95_d[pos]) + '\t' + OT90_offtargets + '\t'
                                       + str(CN_d[pos]) + '\t'
                                       + oligotype + '\t' + FW_oligo + '\t' + RV_oligo + '\t' + res_d[pos] + '\t'
                                       + alias)

                ####################################################
                ##########     continue with guide selection
                ####################################################

                # rescore all other guides based on position relative to new hitzone (distinguish 2nd or 1st hit)
                if hit_type == '2nd':
                    for pos in sgRNAs[gene]:
                        dist = distance_d[pos]
                        if new_zone[0] < dist < new_zone[1]:
                            VBC_score_picking = sgRNAs[gene][pos]
                            VBC_score_picking = VBC_score_picking * penality_3rdsgRNA
                            sgRNAs[gene][pos] = VBC_score_picking
                elif hit_type == '1st':
                    new_zone = [dist-region_window/2,dist+region_window/2]
                    hitzone.append(new_zone)
                    for pos in sgRNAs[gene]:
                        dist = distance_d[pos]
                        if new_zone[0] < dist < new_zone[1]:
                            VBC_score_picking = sgRNAs[gene][pos]
                            VBC_score_picking = VBC_score_picking * penality_2ndsgRNA
                            sgRNAs[gene][pos] = VBC_score_picking
                # eliminate (downscore by factor 0.1) all guides targeting +-2 bp region
                for pos in sgRNAs[gene]:
                    candidate_cutsite = pos_to_cutsite(pos)
                    candidate_orientation = pos.split('(')[1][0]
                    if candidate_orientation == orientation_picked_sgRNA and abs(cutsite_picked_sgRNA - candidate_cutsite) <= 2:
                        VBC_score_picking = sgRNAs[gene][pos]
                        VBC_score_picking = VBC_score_picking * proximity_lessthantwobp_penality
                        sgRNAs[gene][pos] = VBC_score_picking
                # break after 1st guide and rerun until 6 sgRNAs picked (outside loop).
                break
            #print(hitzone)
    #pickle.dump(top6_d, open(Species + 'Okt_top6_d.sav','wb'))

    """
    outfile_all.write('genename' + '\t' + 'position' + '\t' + 'automated selection "top 6"'+ '\t' +
                        'recommended sgRNA-length' + '\t' + 'cutting site amino acids' + '\t' +
                        'distance TSS=0 to stop=1' + '\t' + 'VBC score 0=bad_1=good' + '\t' +
                        'sgRNA_activity 0=bad 1=good' + '\t' + 'inframe prediction inDelphi' + '\t' +
                        'Bioscore 0=bad 1=good' + '\t' + 'very likely off-targets OT>95' + '\t' +
                        'very likely off-target genes' + '\t' + 'likely off-targets OT>90' + '\t' +
                        'likely off-target genes' + '\t' + 'perfect match sites in genome' + '\t' +
                        'GeCKO cloning FW_oligo' + '\t' + 'GeCKO cloning RV_oligo' +
                        'CAUTION - Restriction site in sgRNA')
    """




