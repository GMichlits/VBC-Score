__author__ = 'georg.michlits'

#################################################################################################################
#   USE VBC SCORE python script
#   VERSION 0.0
#   05.04.2020
#   by Georg Michlits
################################################################################################################

Input_File = raw_input('enter input file name (press enter for demo input file): ') or 'demo_input.tsv'

################################################################################################################
#   DATE TABLE CONTAIN PARAMETERS FOR FEATURE COMBINATION
################################################################################################################

print('load coefficients')
# those are the coefficients used to combine features.
Bioscore_coef = [-1.35560617e-02,  7.28768830e-03, -9.10872435e-02, -9.12093203e-02,
  5.20895096e-02,  1.46132540e-02,  5.25554849e-05, -1.70952409e-03,
  4.85115276e-04, -2.91539861e-01, -3.29682194e-03,  6.40238792e-04,
 -4.43685642e-03, -2.84532370e-03, -3.81188449e-03, -2.33663733e-02,
 -1.87065701e-02,  1.58079979e-02,  3.49925940e-03, -2.38165778e-03,
 -1.37843451e-02, -2.70233735e-02,  1.04900183e-02, -4.98751951e-03,
 -2.07810054e-02, -1.44286773e-02, -2.96145725e-02, -3.14324555e-02,
 -5.72708383e-02, -4.34242894e-03, -9.62821065e-03, -4.21773395e-02,
 -4.05963458e-02, -3.85751374e-02,  8.98139887e-03]
Bioscore_intercept = -0.1387530587320981

# those are the coefficients used to combine Doench2016 and tracr_v2 score to give the combined sgRNA score.
combo_sgRNA_coef = [0.90761809, 0.95202036]
combo_sgRNA_intercept = 0.103614766

# those are the coefficients used to combine Frameshift prediction, combined sgRNA score and Bioscore.
VBC_score_coef = [0.04217899, 0.10455383, 0.11984344]
VBC_score_intercept = 0.618668409

# those are score_tables for the AA-identity score, every AA at every position has a predefined score.
AA_score_m2 = {'R':5.846203346,'_':0,'H':4.91773677,'K':-1.285257016,'D':-1.015620198,'E':5.364986841,'S':-4.030242353,
               'T':1.914802324,'N':1.276835696,'Q':2.707695716,'G':-2.45899285,'A':3.287199755,'V':-0.671847454,
               'I':5.762692091,'L':-2.849531339,'M':5.624933291,'F':-5.819300973,'Y':7.851767603,'W':2.23518698,
               'C':4.643814305,'P':-9.723673502}
AA_score_m1 = {'R':4.315638103,'_':0,'H':7.122821897,'K':-6.540897653,'D':1.088937631,'E':-4.873514904,'S':-6.607109094,
               'T':0.092719165,'N':1.145215215,'Q':-8.710620858,'G':-11.78206484,'A':-5.983745207,'V':-1.711103606,
               'I':14.39380474,'L':7.431580596,'M':4.475959969,'F':10.85602466,'Y':19.67026922,'W':13.64411646,
               'C':1.337487544,'P':-3.962570328}
AA_score_0 = {'R':3.101349914,'_':0,'H':10.70349912,'K':-7.298967006,'D':1.777999455,'E':-1.839045145,'S':-3.306735992,
              'T':2.113845781,'N':-4.900001503,'Q':-5.018818155,'G':-13.44639932,'A':-17.11116243,'V':7.369381582,
              'I':15.10297245,'L':7.163197495,'M':14.52657565,'F':-1.542719083,'Y':16.32679676,'W':19.30822843,
              'C':9.501680069,'P':-5.252965723}
AA_score_p1 = {'R':3.180523167,'_':0,'H':-0.752043523,'K':-0.014000872,'D':0.403956122,'E':-1.213473269,
               'S':-13.30322243,'T':-4.781796917,'N':-7.799148233,'Q':-2.015624985,'G':3.73300174,'A':-5.5035969,
               'V':5.503479756,'I':7.461828958,'L':5.711620071,'M':5.451831223,'F':3.508486316,'Y':-0.755052896,
               'W':7.709063885,'C':5.130463965,'P':-18.23390244}
AA_score_p2 = {'R':1.479560774,'_':0,'H':4.046558264,'K':-12.7793265,'D':-3.937236666,'E':-4.697809713,
               'S':-7.260958891,'T':4.47321674,'N':-5.310783271,'Q':-8.943404865,'G':-3.48788691,'A':2.72872104,
               'V':9.232360059,'I':10.15586675,'L':13.27352171,'M':-4.461580754,'F':8.352705947,'Y':6.116877132,
               'W':17.91438276,'C':10.80891866,'P':-8.831802177}

# amino acid letters
all_AA = 'RHKDESTNQGAVILMFYWCP'

########################################################
########################################################
#For future versions
#   The median values are species dependent
#   median values can be used only to substitute missing values or data.
#   Species denpendency can be added optionally.
#   If all features are available they are irrelevant.
Species = 'Hm'
if Species == 'Hm':
    median_phylo = 2.64843
    median_phast = 0.79467
    median_AA = 0.7853857142857142
    median_exlen = 199.0
if Species == 'Ms':
    median_phylo = 2.23205
    median_phast = 0.85024
    median_AA = 0.8273285714285714
    median_exlen = 197.0
if Species == 'dm6':
    median_phylo = 4.7681
    median_phast = 0.735143
    median_AA = 0.6393485714285715
    median_exlen = 730.0
if Species == 'ce11':
    median_phylo = 1.98557
    median_phast = 0.769143
    median_AA = 0.5911257142857143
    median_exlen = 255.0
if Species == 'rn6':
    median_phylo = 1.26919
    median_phast = 0.957238
    median_AA = 0.8368785714285715
    median_exlen = 191.0
if Species == 'xenTro9':
    median_phylo = 0.854857
    median_phast = 0.992143
    median_AA = 0.8967542857142858
    median_exlen = 168.0
median_frameshift = 1-0.2865

#   Values come from Human guides and will be applied to all species for scaling.
#   Scaling is done to put values into the range of 0 - 1.
Bioscore_max = 1.09204
Bioscore_min = -1.02802

#make a Trial sheet for Bioscore calculation
Featuresheet = open(Input_File,'r')
Featuresheet_out = open('output_' + Input_File,'w')

for i,line in enumerate(Featuresheet):
    if i == 0:
        Featuresheet_out.write(line.rstrip('\n') + '\tcombined_sgRNA_Score\tBioscore\tVBC-Score')
    if i > 0:
        Featuresheet_out.write('\n' + line.rstrip('\n'))
        column = line.rstrip('\n').split('\t')
        pos = column[0]
        nt30seq = column[1]
        phylo = float(column[2])
        phast = float(column[3])
        Pfamdom = float(column[4])
        Pfamfam = float(column[5])
        distance = float(column[6])
        mod3 = float(column[7])
        exlen = float(column[8])
        exdist = float(column[9])
        if exdist > 21:
            exdist = 21
        exstop = float(column[10])
        if exstop > 21:
            exstop = 21
        AAcons = float(column[11])
        AAtype_m = column[12]
        AAtype_cut = column[13]
        AAtype_p = column[14]
        Frameshiftrate_inDelphi = float(column[15])
        Doench2016_score = float(column[16])
        tracr_v2_score = float(column[17])
        if 2 <= len(AAtype_m):
            AA_m2 = AAtype_m[len(AAtype_m)::-1][1]
            AA_m2_num = AA_score_m2[AA_m2]
        else:
            AA_m2 = 0
        if 1 <= len(AAtype_m):
            AA_m1 = AAtype_m[len(AAtype_m)::-1][0]

            AA_m1_num = AA_score_m1[AA_m1]
        else:
            AA_m1 = 0
        if AAtype_cut in all_AA and not AAtype_cut == '':
            AA_0_num = AA_score_0[AAtype_cut]
        else:
            AA_0_num = 0
        if len(AAtype_p) >= 1:
            AA_p1 = AAtype_p[0]
            AA_p1_num = AA_score_p1[AA_p1]
        else:
            AA_p1_num = 0
        if len(AAtype_p) >= 2:
            AA_p2 = AAtype_p[1]
            AA_p2_num = AA_score_p2[AA_p2]
        else:
            AA_p2_num = 0
        AAwindow13 = AAtype_m+AAtype_cut+AAtype_p
        AAwindow13_bin = []
        for letter in all_AA:
            if letter in AAwindow13:
                AAwindow13_bin.append(1)
            else:
                AAwindow13_bin.append(0)
        Bioscore_raw = (Bioscore_intercept
                   + Bioscore_coef[0] * phylo + Bioscore_coef[1] * phast + Bioscore_coef[2] * Pfamdom
                   + Bioscore_coef[3] * Pfamfam + Bioscore_coef[4] * distance + Bioscore_coef[5] * mod3
                   + Bioscore_coef[6] * exlen + Bioscore_coef[7] * 21 + Bioscore_coef[8] * 21
                   + Bioscore_coef[9] * AAcons + Bioscore_coef[10] * AA_m2_num + Bioscore_coef[11] * AA_m1_num
                   + Bioscore_coef[12] * AA_0_num + Bioscore_coef[13] * AA_p1_num + Bioscore_coef[14] * AA_p2_num)
        for n,letter in enumerate(all_AA):
            Bioscore_raw += Bioscore_coef[15+n] * AAwindow13_bin[n]
            '\tBioscore_raw\tminmaxBioscore\tscaledBioscore\tdictBioscore'
        # the Bioscore is inverted
        Bioscore_raw = Bioscore_raw * -1
        # the Bioscore is scaled to range from 0 to 1
        minmaxBio = (Bioscore_raw-Bioscore_min)/(Bioscore_max-Bioscore_min)
        if minmaxBio > 0.5:
            scaledBio = (minmaxBio-0.5)/0.5*0.9+0.1
        else:
            scaledBio = (minmaxBio)/0.5*0.1
        Bioscore = scaledBio

        # EXTEND THE INPUT FILE BY combined sgRNA score, write output
        combo_sgRNA = Doench2016_score * combo_sgRNA_coef[0] + tracr_v2_score * combo_sgRNA_coef[1] \
                      + combo_sgRNA_intercept
        Featuresheet_out.write('\t' + str(combo_sgRNA))

        # EXTEND THE INPUT FILE BY BIOSCORE, write output
        Featuresheet_out.write('\t' + str(Bioscore))

        # EXTEND THE INPUT FILE BY VBC score, write output
        VBC_raw = Frameshiftrate_inDelphi * VBC_score_coef[0] + combo_sgRNA * VBC_score_coef[1] \
              + Bioscore * VBC_score_coef[2] + VBC_score_intercept
        VBC = (VBC_raw-0.62226)/(0.88656-0.62226) # scaling to get values of 0 to 1
        Featuresheet_out.write('\t' + str(VBC))
        if i % 100 == 0:
            print(str(i) + ' sgRNAs processed')
print('write output: output_' + Input_File)