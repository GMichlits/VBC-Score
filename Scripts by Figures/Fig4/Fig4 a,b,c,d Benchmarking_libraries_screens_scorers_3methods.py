__author__ = 'georg.michlits'

import matplotlib.pyplot as plt
import numpy as np
import pickle
import os
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
Species = 'Ms'
Analysis_name = Species + '_compare2_mm10_Ms'
# make directory for plots
plot_folder = Analysis_name + '_plots'
if not os.path.isdir(plot_folder):
    os.mkdir(plot_folder)

##################################################################################################################################
#################################################### HUMAN
Screen_name_list = []
Screen_data_d_list = []
Screen_data_all_d_list = []

print('load datasets')
if Species == 'Hm':
    #### first list traing set

    Screen_name_list.append('Munoz NCIH1299')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/Novartis_tiled/NCIH1299_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/Novartis_tiled/NCIH1299_abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('Munoz RKO')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/Novartis_tiled/RKO_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/Novartis_tiled/RKO_abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('Munoz DLD')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/Novartis_tiled/DLD_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/Novartis_tiled/DLD_abs_LFC_d.sav', 'rb')))

    #### then v1 tracr screens
    Screen_name_list.append('Geckov2_HT29')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/GeckoV2_Exp_12_HT29_LARGE_INTESTINE_CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/GeckoV2_Exp_12_HT29_LARGE_INTESTINE_abs_d.sav', 'rb')))
    Screen_name_list.append('Geckov2_NCIH2009')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/GeckoV2_NCIH2009_Exp55CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/GeckoV2_NCIH2009_Exp55abs_LFC_d.sav', 'rb')))

    Screen_name_list.append('Geckov2_K562')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/GeckoV2_Exp_13_K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE_CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/GeckoV2_Exp_13_K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE_abs_d.sav', 'rb')))
    Screen_name_list.append('Geckov2 PC3')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Gecko_v2_Exp_25_PC3_PROSTATE_CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Gecko_v2_Exp_25_PC3_PROSTATE_abs_d.sav', 'rb')))
    Screen_name_list.append('GeckoV2 K562 - Exp49')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/GeckoV2_K562_Exp49CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/GeckoV2_K562_Exp49abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('Hart2015 HTC116')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/TKOv1_HTC116_CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/TKOv1_HTC116_abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('Hart2015 DLD1')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Hart2015_DLD1_Exp4CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Hart2015_DLD1_Exp4abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('Hart2015 RPE1')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Hart2015_RPE1_Exp11CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Hart2015_RPE1_Exp11abs_LFC_d.sav', 'rb')))

    Screen_name_list.append('HartMoffat TKOv3')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/TKOv3_HAP1d18_CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/TKOv3_HAP1d18_abs_d.sav', 'rb')))

    Screen_name_list.append('Wang14_KBM7')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Wang2014_Exp17_KBM7_CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Wang2014_Exp17_KBM7_abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('Wang14_HL60')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Wang2014_Exp16_HL60_CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Wang2014_Exp16_HL60_abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('Wang15_Raji')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Wang15_Raji_Exp2CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Wang15_Raji_Exp2abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('Wang15_KBM7')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Wang2015_KBM7_Exp1CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Wang2015_KBM7_Exp1abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('Wang15_K562')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Wang15_K562_Exp3CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Wang15_K562_Exp3abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('Wang17_MOLM13')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Wang17_MOLM13_Exp136CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Wang17_MOLM13_Exp136abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('Wang17_PL21')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Wang17_PL21_Exp128CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Wang17_PL21_Exp128abs_LFC_d.sav', 'rb')))

    Screen_name_list.append('Avana Karpas')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Avana_Karpas_genCR_Exp285_CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Avana_Karpas_genCR_Exp285_abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('Avana HCC1143')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Avana_HCC1143_genCR_Exp227CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Avana_HCC1143_genCR_Exp227abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('Avana MIAPACA2')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Exp340CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Exp340abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('Avana HCC1806')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Exp233CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Exp233abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('Avana JIMT1')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Exp281CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Exp281abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('Avana SKNAS')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Exp415CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Exp415abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('Avana OVCAR8')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Exp375CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Exp375abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('Avana 639V')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Exp151CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Exp151abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('Avana HMC18')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Exp249CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Exp249abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('Avana JMSU1')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Exp282CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Exp282abs_LFC_d.sav', 'rb')))

#### then v2 tracr screens
    Screen_name_list.append('SangerDepmap_OVR8')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/Sanger_depmap_data+spikein/OVR8_CEG3_abs_LFC.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/Sanger_depmap_data+spikein/OVR8_abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('SangerDepmap_RKO-P1D22')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/Sanger_depmap_data+spikein/RKO-P1D22_CEG3_abs_LFC.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/Sanger_depmap_data+spikein/RKO-P1D22_abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('SangerDepmap_AML2-D21')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/Sanger_depmap_data+spikein/AML2-D21_CEG3_abs_LFC.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/Sanger_depmap_data+spikein/AML2-D21_abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('SangerDepmap_ARH77')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/Sanger_depmap_data+spikein/ARH77_CEG3_abs_LFC.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/Sanger_depmap_data+spikein/ARH77_abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('SangerDepmap_H23')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/Sanger_depmap_data+spikein/H23_CEG3_abs_LFC.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/Sanger_depmap_data+spikein/H23_abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('SangerDepmap_HCT116-P1D22')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/Sanger_depmap_data+spikein/HCT116-P1D22_CEG3_abs_LFC.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/Sanger_depmap_data+spikein/HCT116-P1D22_abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('SangerDepmap_AML2')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/Sanger_depmap_data+spikein/AML2_CEG3_abs_LFC.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/Sanger_depmap_data+spikein/AML2_abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('SangerDepmap_SU10')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/Sanger_depmap_data+spikein/SU10_CEG3_abs_LFC.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/Sanger_depmap_data+spikein/SU10_abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('SangerDepmap_SU8')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/Sanger_depmap_data+spikein/SU8_CEG3_abs_LFC.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/Sanger_depmap_data+spikein/SU8_abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('SangerDepmap_CAL27')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/Sanger_depmap_data+spikein/CAL27_CEG3_abs_LFC.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/Sanger_depmap_data+spikein/CAL27_abs_LFC_d.sav', 'rb')))

    Screen_name_list.append('Yusa HL60 Exp3')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/Yusa_v2/Exp_data_thirdCEG/Exp3third_abs_LFC.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/Yusa_v2/Exp_data/Exp3abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('Yusa HT1080 Exp4')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/Yusa_v2/Exp_data_thirdCEG/Exp4third_abs_LFC.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/Yusa_v2/Exp_data/Exp4abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('Yusa MV411 Exp2')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/Yusa_v2/Exp_data_thirdCEG/Exp2third_abs_LFC.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/Yusa_v2/Exp_data/Exp2abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('Yusa MOLM13 Exp12')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/Yusa_v2/Exp_data_thirdCEG/Exp12third_abs_LFC.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/Yusa_v2/Exp_data/Exp12abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('Yusa HT29')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Exp29CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Exp29abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('Yusa OCIAML2')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Exp24CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Exp24abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('Yusa OCIAML3')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Exp23CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Exp23abs_LFC_d.sav', 'rb')))

    Screen_name_list.append('Munoz DLD1')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Exp18CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Exp18abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('Munoz HT1080')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Exp19CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Exp19abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('Munoz MKN45')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Exp20CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Exp20abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('Munoz RKO')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Exp21CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Exp21abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('Munoz SF268')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Exp22CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Exp22abs_LFC_d.sav', 'rb')))


    Screen_name_list.append('Brunello A375')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Brunello_CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Brunello_abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('IMP KBM7')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/Zuber_v3/KBM7_medm1_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/Zuber_v3/KBM7_abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('IMP RKO')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/Zuber_v3/RKO_medm1_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/Zuber_v3/RKO_abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('IMP MIApaca2')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/Zuber_v3/MIAPACA2_medm1_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/Zuber_v3/MIAPACA2_abs_LFC_d.sav', 'rb')))

##################################################################################################################################
#################################################### MOUSE
else:
    Species = 'Ms'
    Screen_name_list = []
    Screen_data_d_list = []
    Screen_data_all_d_list = []
    Screen_name_list.append('CrUMI mESC Data_2n')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/CrUMI_2n_CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/CrUMI_2n_abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('CrUMI mESC Data_n')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/CrUMI_n_CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/CrUMI_n_abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('IMP mESC Data_2n')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Zub_2n_d18_CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Zub_2n_d18_abs_LFC_d.sav', 'rb')))
    Screen_name_list.append('IMP mESC Data_n')
    Screen_data_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Zub_n_d18_CEG3_abs_LFC_d.sav', 'rb')))
    Screen_data_all_d_list.append(pickle.load(open('../../Gen_data/v3_DICTIONARIES/Zub_n_d18_abs_LFC_d.sav', 'rb')))

#print('load model_name')
#model_T1_AAw5_6prop_abs_LFC = pickle.load(open('/Users/georg.michlits/Desktop/Models/VBC_score/model_Nov_rel.sav','rb'))
# gn30 Ms is gn63 Hm
print('load sgRNA property dictionaries')
D16_AA_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/D16_AA_d.sav', 'rb'))
Hart17_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/Hart17_d.sav', 'rb'))
if Species == 'Hm':
    deepCRISPR_K563_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/deepCRISPR_K562_d.sav', 'rb'))
    CRISPRO_d = pickle.load(open('../../Gen_properties/' + Species + '/property_sav/CRISPRO_d.sav','rb'))
CRISPRater_d = pickle.load(open('../../Gen_properties/' + Species + '/property_sav/CRISPRater_d.sav', 'rb'))
deepCRISPR_nt_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/deepCRISPR_nt_d.sav', 'rb'))
Scorer2_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/Scorer2.0_d.sav', 'rb'))
tracrv2_sub_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/tracrv2_d.sav', 'rb'))
tracrv2_rev_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/tracrv2_3Nov_d.sav', 'rb'))
#
#tracrv2_rev_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/tracrv2_3Nov_trainA_d.sav', 'rb'))
#VBC_score_rev_d_rev_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/tracrv2_3Nov_trainB_d.sav', 'rb'))
#D14_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/tracrv2_3Nov_trainC_d.sav', 'rb'))
#
D14_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/D14_d.sav', 'rb'))
VBC_score_rev_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/Jul14_VBC_d.sav', 'rb'))
Bioscore_rev_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/Jul14_Bioscore_d.sav', 'rb'))
VBC_score_sub_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/Apr8_VBC_score_all_rel_v2_d.sav', 'rb'))
Bioscore_sub_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/Apr8_Bioscore_all_rel_v2_d.sav', 'rb'))
print('load Doench_without AApositional information')
D16_woAA_d = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/D16_woAA_d.sav', 'rb'))
print('load inDelphi')
inDelphi = pickle.load(open('../../Gen_properties/'+Species+'/property_sav/inDelphi_infr_d.sav', 'rb'))

rand_d = {}
combo_sgRNA_d = {}
for pos in D16_AA_d:
    rand_d[pos] = np.random.random()
    combo_sgRNA_d[pos] = (tracrv2_rev_d[pos]*0.95202036 + D16_woAA_d[pos]*0.90761809)

pickle.dump(combo_sgRNA_d, open('../../Gen_properties/' + Species + '/property_sav/combo_sgRNA_d.sav','wb'))

inDelphi_fr_d = {}
for pos in inDelphi:
    inDelphi_fr_d[pos] = 1-inDelphi[pos]

print('load other useful dictionaries')
Gene_d = pickle.load(open('../../Gen_properties/' + Species + '/property_sav/Gene_d.sav', 'rb'))
CEG_d = pickle.load(open('../../Gen_properties/' + Species + '/property_sav/CEG_pos_d.sav', 'rb'))
nonEG_d = pickle.load(open('../../Gen_properties/' + Species + '/property_sav/nonEG_d.sav', 'rb'))
CN_d = pickle.load(open('../../Gen_properties/' + Species + '/property_sav/CN_d.sav', 'rb'))

def det_dAUC_best_worst_depleters_only(Screen_data_d, VBC_score_d, D16_AA_d, data_name, print_plot=False):
    Structure_dict = {}
    Screen_data_d_cur = {}
    for pos in Screen_data_d:
        if pos in VBC_score_d:
            if pos in D16_AA_d:
                Screen_data_d_cur[pos] = Screen_data_d[pos]
                LFC = Screen_data_d[pos]
                gene = Gene_d[pos]
                VBC_score = VBC_score_d[pos]
                D16_score = D16_AA_d[pos]
                if not gene in Structure_dict:
                    Structure_dict[gene] = [[],[],[],[]]
                    ### [0] is pos, [1] doench score, [2] VBC score
                Structure_dict[gene][0].append(pos)
                Structure_dict[gene][1].append(D16_score)
                Structure_dict[gene][2].append(VBC_score)
                Structure_dict[gene][3].append(LFC)
    D16 = {}
    D16['best'] = {}
    D16['worst'] = {}
    VBC = {}
    VBC['best'] = {}
    VBC['worst'] = {}
    Rand = {}
    Rand['best'] = {}
    Rand['worst'] = {}
    Perfect = {}
    Perfect['best'] = {}
    Perfect['worst'] = {}
    n = 0
    guide_per_gene = []
    for gene in Structure_dict:
        n += 1
        triple_list = Structure_dict[gene]
        guide_per_gene.append(len(triple_list[0]))
        D16_best_pos = triple_list[0][triple_list[1].index(np.max(triple_list[1]))]
        Enrico_best_pos = triple_list[0][triple_list[2].index(np.max(triple_list[2]))]
        Perfect_best_pos = triple_list[0][triple_list[3].index(np.min(triple_list[3]))]
        D16_worst_pos = triple_list[0][triple_list[1].index(np.min(triple_list[1]))]
        Enrico_worst_pos = triple_list[0][triple_list[2].index(np.min(triple_list[2]))]
        Perfect_worst_pos = triple_list[0][triple_list[3].index(np.max(triple_list[3]))]
        Rand_worst_pos = triple_list[0][np.random.randint(len(triple_list[0]))]
        Rand_best_pos = triple_list[0][np.random.randint(len(triple_list[0]))]
        D16['best'][D16_best_pos] = 0
        D16['worst'][D16_worst_pos] = 0
        VBC['best'][Enrico_best_pos] = 0
        VBC['worst'][Enrico_worst_pos] = 0
        Rand['best'][Rand_best_pos] = 0
        Rand['worst'][Rand_worst_pos] = 0
        Perfect['best'][Perfect_best_pos] = 0
        Perfect['worst'][Perfect_worst_pos] = 0

    i = 0
    x = []
    EN_best_AC = []
    EN_best = 0
    EN_worst_AC = []
    EN_worst = 0
    D16_best_AC = []
    D16_best = 0
    D16_worst_AC = []
    D16_worst = 0
    Rand_best_AC = []
    Rand_best = 0
    Rand_worst_AC = []
    Rand_worst = 0
    Perfect_best_AC = []
    Perfect_best = 0
    Perfect_worst_AC = []
    Perfect_worst = 0
    for item in sorted(Screen_data_d.items(), key=lambda t: t[1]):
        i += 1
        x.append(i)
        pos = item[0]
        LFC = item[1]
        #print(pos + '\t' + str(LFC))
        if pos in VBC['best']:
            EN_best += 1/n
        if pos in VBC['worst']:
            EN_worst += 1/n
        EN_best_AC.append(EN_best)
        EN_worst_AC.append(EN_worst)
        if pos in D16['best']:
            D16_best += 1/n
        if pos in D16['worst']:
            D16_worst += 1/n
        D16_best_AC.append(D16_best)
        D16_worst_AC.append(D16_worst)
        if pos in Rand['best']:
            Rand_best += 1/n
        if pos in Rand['worst']:
            Rand_worst += 1/n
        Rand_best_AC.append(Rand_best)
        Rand_worst_AC.append(Rand_worst)
        if pos in Perfect['best']:
            Perfect_best += 1/n
        if pos in Perfect['worst']:
            Perfect_worst += 1/n
        Perfect_best_AC.append(Perfect_best)
        Perfect_worst_AC.append(Perfect_worst)

    if print_plot:
        plt.title('test on' + data_name)
        plt.ylabel('AC_discovery')
        plt.xlabel('guide_ranking')
        plt.scatter(x, EN_best_AC, c='#FF00FF', s=1)
        plt.scatter(x, EN_worst_AC, c='#FF00FF', s=1)
        plt.scatter(x, D16_best_AC, c='#008000', s=1)
        plt.scatter(x, D16_worst_AC, c='#008000', s=1)
        plt.scatter(x, Rand_best_AC, c='black', s=1)
        plt.scatter(x, Rand_worst_AC, c='black', s=1)
        plt.scatter(x, Perfect_best_AC, c='grey', s=1)
        plt.scatter(x, Perfect_worst_AC, c='grey', s=1)
        plt.savefig(Analysis_name + '_' + data_name +'_AUC.png')
        plt.close()
        plt.close()
        plt.close()
        plt.close()
        plt.close()
        plt.close()


    dAUC_VBC = (sum(EN_best_AC)-sum(EN_worst_AC))/i
    dAUC_D16 = (sum(D16_best_AC)-sum(D16_worst_AC))/i
    dAUC_rand = (sum(Rand_best_AC)-sum(Rand_worst_AC))/i
    dAUC_Perfect = (sum(Perfect_best_AC)-sum(Perfect_worst_AC))/i
    return (dAUC_VBC,dAUC_D16,dAUC_rand,dAUC_Perfect)
def det_dAUC_best_worst_depleters_only_20percent(Screen_data_d, VBC_score_d, D16_AA_d, data_name, print_plot=False):
    Structure_dict = {}
    Screen_data_d_cur = {}
    for pos in Screen_data_d:
        if pos in VBC_score_d:
            if pos in D16_AA_d:
                Screen_data_d_cur[pos] = Screen_data_d[pos]
                LFC = Screen_data_d[pos]
                gene = Gene_d[pos]
                VBC_score = VBC_score_d[pos]
                D16_score = D16_AA_d[pos]
                if not gene in Structure_dict:
                    Structure_dict[gene] = [[],[],[],[]]
                    ### [0] is pos, [1] doench score, [2] VBC score [3] is LFC
                Structure_dict[gene][0].append(pos)
                Structure_dict[gene][1].append(D16_score)
                Structure_dict[gene][2].append(VBC_score)
                Structure_dict[gene][3].append(LFC)
    D16 = {}
    D16['best'] = {}
    D16['worst'] = {}
    VBC = {}
    VBC['best'] = {}
    VBC['worst'] = {}
    Rand = {}
    Rand['best'] = {}
    Rand['worst'] = {}
    Perfect = {}
    Perfect['best'] = {}
    Perfect['worst'] = {}
    n = 0
    guide_per_gene = []
    for gene in Structure_dict:
        #initial_length_gene = len(Structure_dict[gene][0])
        #termination_length = int(initial_length_gene*0.6)
        #60% because initially it shoudl detect the best and worst sgRNA per gene
        #and for e.g. 5 sgRNAs 20% best and 20% worst would be 1 sgRNA selected each
        # for libraries with 10sgRNAs 20% would be top 2 bottom 2
        # for tiled dataset with 100 sgRNAs it would selecet 20 top and 20 bad sgRNAs...
        if len(Structure_dict[gene][0]) >= 3:
            run_n = 0
            d16_pos_list = Structure_dict[gene][0].copy()
            d16_d16_list = Structure_dict[gene][1].copy()
            VBC_pos_list = Structure_dict[gene][0].copy()
            VBC_VBC_list = Structure_dict[gene][2].copy()
            rand_pos_list = Structure_dict[gene][0].copy()
            #for random no lookup-list required
            perf_pos_list = Structure_dict[gene][0].copy()
            perf_LFC_list = Structure_dict[gene][3].copy()
            ########################
            # make 20% d16 list
            initial_len = len(d16_pos_list)
            termin_len = initial_len * 0.6
            guide_per_gene.append(initial_len)
            #extract top/bottom D16
            while len(VBC_pos_list) > termin_len:
                n += 1  #this n is important!! - it count total number of sgRNAs that are in best/worst ranking
                        # and is required to calculate cumulative fraction e.g. 500 sgRNA out of 1000 best sgRNAs
                        # found gives 500/n cumulative fraction (y-axis)
                max_index_VBC = VBC_VBC_list.index(np.max(VBC_VBC_list))
                VBC_best_pos = VBC_pos_list[max_index_VBC]
                VBC['best'][VBC_best_pos] = 0
                VBC_pos_list.pop(max_index_VBC)
                VBC_VBC_list.pop(max_index_VBC)

                min_index_VBC = VBC_VBC_list.index(np.min(VBC_VBC_list))
                VBC_worst_pos = VBC_pos_list[min_index_VBC]
                VBC['worst'][VBC_worst_pos] = 0
                VBC_pos_list.pop(min_index_VBC)
                VBC_VBC_list.pop(min_index_VBC)
            while len(d16_pos_list) > termin_len:
                max_index_D16 = d16_d16_list.index(np.max(d16_d16_list))
                d16_best_pos = d16_pos_list[max_index_D16]
                D16['best'][d16_best_pos] = 0
                d16_pos_list.pop(max_index_D16)
                d16_d16_list.pop(max_index_D16)

                min_index_D16 = d16_d16_list.index(np.min(d16_d16_list))
                d16_worst_pos = d16_pos_list[min_index_D16]
                D16['worst'][d16_worst_pos] = 0
                d16_pos_list.pop(min_index_D16)
                d16_d16_list.pop(min_index_D16)
            #extract top/bottom VBC
            #extract top/bottom sgRNA by random selection random
            while len(rand_pos_list) > termin_len:
                rand_index_VBC = np.random.randint(len(rand_pos_list))
                rand_best_pos = rand_pos_list[rand_index_VBC]
                Rand['best'][rand_best_pos] = 0
                rand_pos_list.pop(rand_index_VBC)

                rand_index_VBC = np.random.randint(len(rand_pos_list))
                rand_worst_pos = rand_pos_list[rand_index_VBC]
                Rand['worst'][rand_worst_pos] = 0
                rand_pos_list.pop(rand_index_VBC)
            # extract top/bottom sgRNA by LFC (hypothetical perfect)
            while len(perf_pos_list) > termin_len:
                min_index_perf = perf_LFC_list.index(np.min(perf_LFC_list))
                perf_best_pos = perf_pos_list[min_index_perf]
                Perfect['best'][perf_best_pos] = 0
                perf_pos_list.pop(min_index_perf)
                perf_LFC_list.pop(min_index_perf)

                max_index_perf = perf_LFC_list.index(np.max(perf_LFC_list))
                perf_worst_pos = perf_pos_list[max_index_perf]
                Perfect['worst'][perf_worst_pos] = 0
                perf_pos_list.pop(max_index_perf)
                perf_LFC_list.pop(max_index_perf)

    i = 0
    x = []
    EN_best_AC = []
    EN_best = 0
    EN_worst_AC = []
    EN_worst = 0
    D16_best_AC = []
    D16_best = 0
    D16_worst_AC = []
    D16_worst = 0
    Rand_best_AC = []
    Rand_best = 0
    Rand_worst_AC = []
    Rand_worst = 0
    Perfect_best_AC = []
    Perfect_best = 0
    Perfect_worst_AC = []
    Perfect_worst = 0
    for item in sorted(Screen_data_d.items(), key=lambda t: t[1]):
        i += 1
        x.append(i)
        pos = item[0]
        LFC = item[1]
        #print(pos + '\t' + str(LFC))
        if pos in VBC['best']:
            EN_best += 1/n
        if pos in VBC['worst']:
            EN_worst += 1/n
        EN_best_AC.append(EN_best)
        EN_worst_AC.append(EN_worst)
        if pos in D16['best']:
            D16_best += 1/n
        if pos in D16['worst']:
            D16_worst += 1/n
        D16_best_AC.append(D16_best)
        D16_worst_AC.append(D16_worst)
        if pos in Rand['best']:
            Rand_best += 1/n
        if pos in Rand['worst']:
            Rand_worst += 1/n
        Rand_best_AC.append(Rand_best)
        Rand_worst_AC.append(Rand_worst)
        if pos in Perfect['best']:
            Perfect_best += 1/n
        if pos in Perfect['worst']:
            Perfect_worst += 1/n
        Perfect_best_AC.append(Perfect_best)
        Perfect_worst_AC.append(Perfect_worst)

    if print_plot:
        plt.title('test on' + data_name)
        plt.ylabel('AC_discovery')
        plt.xlabel('guide_ranking')
        plt.scatter(x, EN_best_AC, c='#FF00FF', s=1)
        plt.scatter(x, EN_worst_AC, c='#FF00FF', s=1)
        plt.scatter(x, D16_best_AC, c='#008000', s=1)
        plt.scatter(x, D16_worst_AC, c='#008000', s=1)
        plt.scatter(x, Rand_best_AC, c='black', s=1)
        plt.scatter(x, Rand_worst_AC, c='black', s=1)
        plt.scatter(x, Perfect_best_AC, c='grey', s=1)
        plt.scatter(x, Perfect_worst_AC, c='grey', s=1)
        plt.savefig(Analysis_name + '_' + data_name + '_'+'AUC.png')
        plt.close()
        plt.close()
        plt.close()
        plt.close()
        plt.close()
        plt.close()


    dAUC_VBC = (sum(EN_best_AC)-sum(EN_worst_AC))/i
    dAUC_D16 = (sum(D16_best_AC)-sum(D16_worst_AC))/i
    dAUC_rand = (sum(Rand_best_AC)-sum(Rand_worst_AC))/i
    dAUC_Perfect = (sum(Perfect_best_AC)-sum(Perfect_worst_AC))/i
    return (dAUC_VBC,dAUC_D16,dAUC_rand,dAUC_Perfect)
def det_dAUC_best_worst_all(Screen_data_all_d, VBC_score_d, D16_AA_d, data_name, model_name, print_plot=False):
    Structure_dict = {}
    Screen_data_d_cur = {}
    for pos in Screen_data_all_d:
        if pos in VBC_score_d:
            if pos in D16_AA_d:
                Screen_data_d_cur[pos] = Screen_data_d[pos]
                gene = Gene_d[pos]
                VBC_score = VBC_score_d[pos]
                D16_score = D16_AA_d[pos]
                if not gene in Structure_dict:
                    Structure_dict[gene] = [[],[],[]]
                    ### [0] is pos, [1] doench score, [2] VBC score
                Structure_dict[gene][0].append(pos)
                Structure_dict[gene][1].append(D16_score)
                Structure_dict[gene][2].append(VBC_score)
    D16 = {}
    D16['best'] = {}
    D16['worst'] = {}
    VBC = {}
    VBC['best'] = {}
    VBC['worst'] = {}
    Rand = {}
    Rand['best'] = {}
    Rand['worst'] = {}
    n = 0
    guide_per_gene = []
    for gene in Structure_dict:
        n += 1
        triple_list = Structure_dict[gene]
        guide_per_gene.append(len(triple_list[0]))
        D16_best_pos = triple_list[0][triple_list[1].index(np.max(triple_list[1]))]
        Enrico_best_pos = triple_list[0][triple_list[2].index(np.max(triple_list[2]))]
        D16_worst_pos = triple_list[0][triple_list[1].index(np.min(triple_list[1]))]
        Enrico_worst_pos = triple_list[0][triple_list[2].index(np.min(triple_list[2]))]
        Rand_worst_pos = triple_list[0][np.random.randint(len(triple_list[0]))]
        Rand_best_pos = triple_list[0][np.random.randint(len(triple_list[0]))]
        D16['best'][D16_best_pos] = 0
        D16['worst'][D16_worst_pos] = 0
        VBC['best'][Enrico_best_pos] = 0
        VBC['worst'][Enrico_worst_pos] = 0
        Rand['best'][Rand_best_pos] = 0
        Rand['worst'][Rand_worst_pos] = 0

    i = 0
    x = []
    EN_best_AC = []
    EN_best = 0
    EN_worst_AC = []
    EN_worst = 0
    D16_best_AC = []
    D16_best = 0
    D16_worst_AC = []
    D16_worst = 0
    Rand_best_AC = []
    Rand_best = 0
    Rand_worst_AC = []
    Rand_worst = 0
    for item in sorted(Screen_data_all_d.items(), key=lambda t: t[1]):
        i += 1
        x.append(i)
        pos = item[0]
        LFC = item[1]
        #print(pos + '\t' + str(LFC))
        if pos in VBC['best']:
            EN_best += 1/n
        if pos in VBC['worst']:
            EN_worst += 1/n
        EN_best_AC.append(EN_best)
        EN_worst_AC.append(EN_worst)
        if pos in D16['best']:
            D16_best += 1/n
        if pos in D16['worst']:
            D16_worst += 1/n
        D16_best_AC.append(D16_best)
        D16_worst_AC.append(D16_worst)
        if pos in Rand['best']:
            Rand_best += 1/n
        if pos in Rand['worst']:
            Rand_worst += 1/n
        Rand_best_AC.append(Rand_best)
        Rand_worst_AC.append(Rand_worst)


    if print_plot:
        plt.title('Train with' + model_name + 'test on' + data_name)
        plt.ylabel('AC_discovery')
        plt.xlabel('guide_ranking')
        plt.scatter(x, EN_best_AC, c='#FF00FF', s=1)
        plt.scatter(x, EN_worst_AC, c='#FF00FF', s=1)
        plt.scatter(x, D16_best_AC, c='#008000', s=1)
        plt.scatter(x, D16_worst_AC, c='#008000', s=1)
        plt.scatter(x, Rand_best_AC, c='black', s=1)
        plt.scatter(x, Rand_worst_AC, c='black', s=1)
        plt.savefig(data_name + '_' + model_name +'AUC.png')
        plt.close()
        plt.close()
        plt.close()
        plt.close()
        plt.close()
        plt.close()


    dAUC_VBC = (sum(EN_best_AC)-sum(EN_worst_AC))/i
    dAUC_D16 = (sum(D16_best_AC)-sum(D16_worst_AC))/i
    dAUC_rand = (sum(Rand_best_AC)-sum(Rand_worst_AC))/i
    return (dAUC_VBC,dAUC_D16,dAUC_rand)
def det_dAUC(Screen_data_all_d, CEG_d, nonEG_d):

    sum_total = 0
    sum_CEG = 0
    sum_nonEG = 0
    for pos_item in sorted(Screen_data_all_d.items(), key=lambda t: t[1]):
        pos = pos_item[0]
        LFC = pos_item[1]
        sum_total += 1
        if pos in CEG_d:
            sum_CEG += 1
        if pos in nonEG_d:
            sum_nonEG += 1
    AC_CEG = 0
    AC_nonEG = 0
    accum_CEG = 0
    accum_nonEG = 0
    x = []
    y_dia = []
    y_CEG = []
    y_nonEG = []
    i = 0
    for pos_item in sorted(Screen_data_all_d.items(), key=lambda t: t[1]):
        i += 1
        pos = pos_item[0]
        LFC = pos_item[1]
        x.append(i)
        y_dia.append(i/sum_total)
        if pos in CEG_d:
            AC_CEG += 1/sum_CEG
        if pos in nonEG_d:
            AC_nonEG += 1/sum_nonEG
        y_CEG.append(AC_CEG)
        accum_CEG += AC_CEG
        y_nonEG.append(AC_nonEG)
        accum_nonEG += AC_nonEG
    total_area = len(x)*1
    AUC_CEG = round(accum_CEG/total_area, 4)
    AUC_nonEG = round(accum_nonEG/total_area, 4)
    dAUC = round(AUC_CEG-AUC_nonEG, 4)
    return(dAUC,AUC_CEG,AUC_nonEG,sum_total,sum_CEG,sum_nonEG,x,y_dia,y_CEG,y_nonEG)
def det_CEG(Screen_data_all_d,CEG_d):
    CEG_list = []
    for pos in Screen_data_all_d:
        if pos in CEG_d:
            CEG_list.append(Screen_data_all_d[pos])
    CEG_depl = np.median(CEG_list)
    return(CEG_depl)
def Pearson_plot_score_dependency_lfc(screen_d, score_d, Copynum_d=CN_d, make_plot=False, savepltfolder = 'folder',
                                      title='', x_axis='score', y_axis='screen LFC', set_enrichers0=True):
    x = []
    y = []
    for pos in screen_d:
        if Copynum_d[pos] == 1:
            if pos in score_d:
                LFC = screen_d[pos]
                score = score_d[pos]
                if set_enrichers0:
                    if LFC > 0:
                        LFC = 0
                x.append(score)
                y.append(LFC)
    PearsCoeff, pears_pval = stats.pearsonr(x,y)
    if make_plot:
        plt.rc('xtick',labelsize=20)
        plt.rc('ytick',labelsize=20)
        plt.figure(figsize=(20,15))
        plt.title(title, size=30)
        plt.scatter(x, y, c='k', s=50,alpha=0.1)
        gradient, intercept, r_value, p_value, std_err = stats.linregress(x,y)
        plt.xlabel(x_axis, size=30)
        plt.ylabel(y_axis, size=30)
        mn = np.min(x)
        mx = np.max(x)
        x1 = np.linspace(mn, mx, 500)
        y1 = gradient*x1+intercept
        plt.plot(x1, y1, c='r', linewidth=2)
        plt.savefig(savepltfolder + '/' + title + '.png')
        plt.close()
    return PearsCoeff, pears_pval

def Pearson_plot_score_dependency_linear(screen_d, score_d, Copynum_d=CN_d, make_plot=False, savepltfolder = 'folder',
                                         title='', x_axis ='score', y_axis ='screen depletion efficiency',
                                         set_enrichers0=True):
    x = []
    y = []
    for pos in screen_d:
        if Copynum_d[pos] == 1:
            if pos in score_d:
                depletion_efficiency = 1-(2**screen_d[pos])
                if set_enrichers0:
                    if depletion_efficiency<0:
                        depletion_efficiency = 0
                score = score_d[pos]
                x.append(score)
                y.append(depletion_efficiency)
    PearsCoeff, pears_pval = stats.pearsonr(x,y)
    if make_plot:
        plt.rc('xtick',labelsize=20)
        plt.rc('ytick',labelsize=20)
        plt.figure(figsize=(20,15))
        plt.title(title, size=30)
        plt.scatter(x,y,c='k',s=50,alpha=0.1)
        gradient, intercept, r_value, p_value, std_err = stats.linregress(x,y)
        plt.xlabel(x_axis, size=30)
        plt.ylabel(y_axis, size=30)
        mn = np.min(x)
        mx = np.max(x)
        x1 = np.linspace(mn, mx, 500)
        y1 = gradient*x1+intercept
        plt.plot(x1, y1, c='r', linewidth=2)
        plt.savefig(savepltfolder + '/' + title + '.png')
        plt.close()

    return PearsCoeff, pears_pval
def split_screen_by_score(Screen_d, score_d, percentage=50, Gene_namefrompos_d=Gene_d, reference_essential_genes_d=CEG_d,
                 reference_non_essentials_d=nonEG_d, Copynumber_dict=CN_d):
    Structure_dict = {}
    for pos in Screen_d:
        if pos in Gene_namefrompos_d and pos in score_d and Copynumber_dict[pos] == 1:
            Gene = Gene_namefrompos_d[pos]
            if not Gene in Structure_dict:
                Structure_dict[Gene] = [[],[],[]] #pos, score, LFC
            Structure_dict[Gene][0].append(pos)
            Structure_dict[Gene][1].append(score_d[pos])
            Structure_dict[Gene][2].append(Screen_d[pos])
    CEG_VBC = {}
    CEG_VBC['best'] = {}
    CEG_VBC['worst'] = {}
    CEG_Rand = {}
    CEG_Rand['best'] = {}
    CEG_Rand['worst'] = {}
    CEG_Perfect = {}
    CEG_Perfect['best'] = {}
    CEG_Perfect['worst'] = {}

    nEG_VBC = {}
    nEG_VBC['best'] = {}
    nEG_VBC['worst'] = {}
    nEG_Rand = {}
    nEG_Rand['best'] = {}
    nEG_Rand['worst'] = {}
    nEG_Perfect = {}
    nEG_Perfect['best'] = {}
    nEG_Perfect['worst'] = {}
    n = 0
    guide_per_gene = []
    # for CEG_selected:
    for gene in Structure_dict:
        first_pos = Structure_dict[gene][0][0]
        if first_pos in reference_essential_genes_d:
            if len(Structure_dict[gene][0]) >= 3:
                run_n = 0
                VBC_pos_list = Structure_dict[gene][0].copy()
                VBC_VBC_list = Structure_dict[gene][1].copy()
                rand_pos_list = Structure_dict[gene][0].copy()
                #rand_rand_list not required
                perf_pos_list = Structure_dict[gene][0].copy()
                perf_LFC_list = Structure_dict[gene][2].copy()
                ########################
                # make 20% d16 list
                initial_len = len(VBC_pos_list)
                termin_len = initial_len * 1-2*(percentage/100)
                guide_per_gene.append(initial_len)
                #extract top/bottom VBC
                while len(VBC_pos_list) > max([termin_len, 1.1]):
                    n += 1  #this n is important!! - it count total number of sgRNAs that are in best/worst ranking
                            # and is required to calculate cumulative fraction e.g. 500 sgRNA out of 1000 best sgRNAs
                            # found gives 500/n cumulative fraction (y-axis)
                    max_index_VBC = VBC_VBC_list.index(np.max(VBC_VBC_list))
                    VBC_best_pos = VBC_pos_list[max_index_VBC]
                    CEG_VBC['best'][VBC_best_pos] = 0
                    VBC_pos_list.pop(max_index_VBC)
                    VBC_VBC_list.pop(max_index_VBC)

                    min_index_VBC = VBC_VBC_list.index(np.min(VBC_VBC_list))
                    VBC_worst_pos = VBC_pos_list[min_index_VBC]
                    CEG_VBC['worst'][VBC_worst_pos] = 0
                    VBC_pos_list.pop(min_index_VBC)
                    VBC_VBC_list.pop(min_index_VBC)
                #extract top/bottom sgRNA by random selection random
                while len(rand_pos_list) > max([termin_len, 1.1]):
                    rand_index_VBC = np.random.randint(len(rand_pos_list))
                    rand_best_pos = rand_pos_list[rand_index_VBC]
                    CEG_Rand['best'][rand_best_pos] = 0
                    rand_pos_list.pop(rand_index_VBC)

                    rand_index_VBC = np.random.randint(len(rand_pos_list))
                    rand_worst_pos = rand_pos_list[rand_index_VBC]
                    CEG_Rand['worst'][rand_worst_pos] = 0
                    rand_pos_list.pop(rand_index_VBC)
                # extract top/bottom sgRNA by LFC (hypothetical perfect)
                while len(perf_pos_list) > max([termin_len, 1.1]):
                    min_index_perf = perf_LFC_list.index(np.min(perf_LFC_list))
                    perf_best_pos = perf_pos_list[min_index_perf]
                    CEG_Perfect['best'][perf_best_pos] = 0
                    perf_pos_list.pop(min_index_perf)
                    perf_LFC_list.pop(min_index_perf)

                    max_index_perf = perf_LFC_list.index(np.max(perf_LFC_list))
                    perf_worst_pos = perf_pos_list[max_index_perf]
                    CEG_Perfect['worst'][perf_worst_pos] = 0
                    perf_pos_list.pop(max_index_perf)
                    perf_LFC_list.pop(max_index_perf)
        elif first_pos in reference_non_essentials_d:
            if len(Structure_dict[gene][0]) >= 3:
                run_n = 0
                VBC_pos_list = Structure_dict[gene][0].copy()
                VBC_VBC_list = Structure_dict[gene][1].copy()
                rand_pos_list = Structure_dict[gene][0].copy()
                #rand_rand_list not required
                perf_pos_list = Structure_dict[gene][0].copy()
                perf_LFC_list = Structure_dict[gene][2].copy()
                ########################
                # make 20% d16 list
                initial_len = len(VBC_pos_list)
                termin_len = initial_len * 1-2*(percentage/100)
                guide_per_gene.append(initial_len)
                #extract top/bottom VBC
                while len(VBC_pos_list) > max([termin_len, 1.1]):
                    n += 1  #this n is important!! - it count total number of sgRNAs that are in best/worst ranking
                            # and is required to calculate cumulative fraction e.g. 500 sgRNA out of 1000 best sgRNAs
                            # found gives 500/n cumulative fraction (y-axis)
                    max_index_VBC = VBC_VBC_list.index(np.max(VBC_VBC_list))
                    VBC_best_pos = VBC_pos_list[max_index_VBC]
                    nEG_VBC['best'][VBC_best_pos] = 0
                    VBC_pos_list.pop(max_index_VBC)
                    VBC_VBC_list.pop(max_index_VBC)

                    min_index_VBC = VBC_VBC_list.index(np.min(VBC_VBC_list))
                    VBC_worst_pos = VBC_pos_list[min_index_VBC]
                    nEG_VBC['worst'][VBC_worst_pos] = 0
                    VBC_pos_list.pop(min_index_VBC)
                    VBC_VBC_list.pop(min_index_VBC)
                #extract top/bottom sgRNA by random selection random
                while len(rand_pos_list) > max([termin_len, 1.1]):
                    rand_index_VBC = np.random.randint(len(rand_pos_list))
                    rand_best_pos = rand_pos_list[rand_index_VBC]
                    nEG_Rand['best'][rand_best_pos] = 0
                    rand_pos_list.pop(rand_index_VBC)

                    rand_index_VBC = np.random.randint(len(rand_pos_list))
                    rand_worst_pos = rand_pos_list[rand_index_VBC]
                    nEG_Rand['worst'][rand_worst_pos] = 0
                    rand_pos_list.pop(rand_index_VBC)
                # extract top/bottom sgRNA by LFC (hypothetical perfect)
                while len(perf_pos_list) > max([termin_len, 1.1]):
                    min_index_perf = perf_LFC_list.index(np.min(perf_LFC_list))
                    perf_best_pos = perf_pos_list[min_index_perf]
                    nEG_Perfect['best'][perf_best_pos] = 0
                    perf_pos_list.pop(min_index_perf)
                    perf_LFC_list.pop(min_index_perf)

                    max_index_perf = perf_LFC_list.index(np.max(perf_LFC_list))
                    perf_worst_pos = perf_pos_list[max_index_perf]
                    nEG_Perfect['worst'][perf_worst_pos] = 0
                    perf_pos_list.pop(max_index_perf)
                    perf_LFC_list.pop(max_index_perf)
    return CEG_VBC, CEG_Rand, CEG_Perfect, nEG_VBC, nEG_Rand, nEG_Perfect

if Species == 'Hm':
    print('data\tmodel\tperfect\tVBC\tD16AA\trand\tBio\tinDeplphi\tCRIPRO\tCRISPRater\tD16woAA\trand\td_AUCess\tCEGdepl\tHart17sc\tD14score\tScorer2.0\ttracrv2\tdeepCRISPR_nt\tdeepCRISPR_K562')
else:
    print('data\tmodel\tperfect\tVBC\tD16AA\trand\tBio\tinDeplphi\tCRISPRO\tCRISPRater\tD16woAA\trand\td_AUCess\tCEGdepl\tHart17sc\tD14score\tScorer2.0\ttracrv2\tdeepCRISPR_nt\tdeepCRISPR_K562')
outfile_dAUC_best_worst = open(Analysis_name + '.txt','w')
outfile_dAUC_best_worst.write('data\tperfect\trand\tD14\tScorer2.0\tHart17sc\tD16woAA\tD16wAA\ttracr_v2_sub\tinDelphi'+
            '\tBioscore_sub\tVBC_score_sub\ttracrv2_rev\tcombo_sgRNA\tinDelphi\tBioscore_rev\tVBC_score_rev'+
            '\tCRISPRater\tCRISPRO\tdeepCRISPR_nt\tdeepCRISPR_K562\td_AUCess\tCEGdepl')
pearson_r_lfc_file = open(Analysis_name + 'pearson_lfc.txt', 'w')
pearson_r_lin_file = open(Analysis_name + 'pearson_lin.txt', 'w')
ddAUC_output = open(Analysis_name + 'ddAUC_ROC_trad.txt', 'w')
rel_to_perfect_ddAUC_output = open(Analysis_name + 'rel_to_perf_ddAUC_ROC_trad.txt','w')

if Species == 'Hm':
    score_list = [rand_d, D14_d, Scorer2_d, Hart17_d, D16_woAA_d, D16_AA_d, tracrv2_sub_d, inDelphi_fr_d, Bioscore_sub_d,
              VBC_score_sub_d, tracrv2_rev_d, combo_sgRNA_d, inDelphi_fr_d, Bioscore_rev_d, VBC_score_rev_d, CRISPRater_d, CRISPRO_d,
              deepCRISPR_nt_d, deepCRISPR_K563_d]
    score_list_names = ['rand', 'Doench 2014', 'Scorer2.0', 'Hart2017', 'Doench 2016,noAA', 'Doench 2016,+AA', 'tracrv2_sub',
                        'inDelphi_fr', 'Bioscore_sub', 'VBC_score_sub', 'tracrv2_rev','combo_sgRNA','inDelphi_fr', 'Bioscore_rev',
                        'VBC_score_rev', 'CRISPRater', 'CRISPRO', 'deepCRISPR_nt', 'eepCRISPR_K562']
if Species == 'Ms':
    score_list = [rand_d, D14_d, Scorer2_d, Hart17_d, D16_woAA_d, D16_AA_d, tracrv2_sub_d, inDelphi_fr_d, Bioscore_sub_d,
              VBC_score_sub_d, tracrv2_rev_d, combo_sgRNA_d, inDelphi_fr_d, Bioscore_rev_d, VBC_score_rev_d, CRISPRater_d, CRISPRater_d,
              deepCRISPR_nt_d, deepCRISPR_nt_d]
    score_list_names = ['rand', 'Doench 2014', 'Scorer2.0', 'Hart2017', 'Doench 2016,noAA', 'Doench 2016,+AA', 'tracrv2_sub',
                        'inDelphi_fr', 'Bioscore_sub', 'VBC_score_sub', 'tracrv2_rev','combo_sgRNA','inDelphi_fr', 'Bioscore_rev',
                        'VBC_score_rev', 'CRISPRater', 'CRISPRater', 'deepCRISPR_nt', 'deepCRISPR_nt']

pearson_r_lfc_file.write('data')
pearson_r_lin_file.write('data')
ddAUC_output.write('data')
rel_to_perfect_ddAUC_output.write('data')
for score in score_list_names:
    pearson_r_lfc_file.write('\t'+score)
    pearson_r_lin_file.write('\t'+score)
    ddAUC_output.write('\t' + score)
    rel_to_perfect_ddAUC_output.write('\t' + score)

x = []
y_Bio = []
y_inD = []
y_D16 = []
y_D16wo = []
y_VBC = []
y_TR2 = []
y_rand = []
c_Bio = []
c_inD = []
c_D16 = []
c_D16wo = []
c_VBC = []
c_TR2 = []
c_rand = []
############################# For EVERY SCREEN LOADED.....
for i in range(len(Screen_name_list)):
    data_name = Screen_name_list[i]
    Screen_data_d = Screen_data_d_list[i]
    Screen_data_all_d = Screen_data_all_d_list[i]
    print('processing ' + data_name + ' ' + str(i+1) + ' out of ' + str(len(Screen_name_list)))

############################# Calculate dAUC (best versus worst 20% of sgRNAs by meassure of the different scores:

    dAUC_VBC_rev,dAUC_D16wAA,dAUC_rand1,dAUC_Perfect = det_dAUC_best_worst_depleters_only_20percent(Screen_data_d, VBC_score_rev_d, D16_AA_d, data_name, print_plot=False)
    dAUC_Bio_rev,dAUC_inDelphi,dAUC_rand2,dAUC_Perfect = det_dAUC_best_worst_depleters_only_20percent(Screen_data_d, Bioscore_rev_d, inDelphi_fr_d, data_name)
    dAUC_VBC_sub,dAUC_Bio_sub,dAUC_rand2,dAUC_Perfect = det_dAUC_best_worst_depleters_only_20percent(Screen_data_d, VBC_score_sub_d, Bioscore_sub_d, data_name)
    dAUC_tracrv2_rev,dAUC_tracrv2_sub,dAUC_rand2,dAUC_Perfect = det_dAUC_best_worst_depleters_only_20percent(Screen_data_d, tracrv2_rev_d, tracrv2_sub_d, data_name)
    dAUC_combo_sgRNA,dAUC_tracrv2_sub,dAUC_rand2,dAUC_Perfect = det_dAUC_best_worst_depleters_only_20percent(Screen_data_d, combo_sgRNA_d, tracrv2_sub_d, data_name)
    dAUC_deepCRISPR_nt,dAUC_CRISPRater,dAUC_rand2,dAUC_Perfect = det_dAUC_best_worst_depleters_only_20percent(Screen_data_d, deepCRISPR_nt_d, CRISPRater_d, data_name)
    if Species == 'Hm':
        dAUC_deepCRISPR_nt,dAUC_deepCRISPR_K562,dAUC_rand2,dAUC_Perfect = det_dAUC_best_worst_depleters_only_20percent(Screen_data_d, deepCRISPR_nt_d, deepCRISPR_K563_d, data_name)
        dAUC_CRISPRater,dAUC_CRISPRO,dAUC_rand2,dAUC_Perfect = det_dAUC_best_worst_depleters_only_20percent(Screen_data_d, CRISPRater_d, CRISPRO_d, data_name)
    dAUC_Hart17,dAUC_D16woAA,dAUC_rand3,dAUC_Perfect = det_dAUC_best_worst_depleters_only_20percent(Screen_data_d, Hart17_d, D16_woAA_d, data_name)
    dAUC_D14,dAUC_Scorer2,dAUC_rand3,dAUC_Perfect = det_dAUC_best_worst_depleters_only_20percent(Screen_data_d, D14_d, Scorer2_d, data_name)
    dAUC_av_rand = np.mean([dAUC_rand1,dAUC_rand2,dAUC_rand3])
    dAUC,AUC_CEG,AUC_nonEG,sum_total,sum_CEG,sum_nonEG,x,y_dia,y_CEG,y_nonEG = det_dAUC(Screen_data_all_d, CEG_d, nonEG_d)
    CEG_depl = det_CEG(Screen_data_all_d, CEG_d)

    if Species == 'Hm':
        line = (data_name + '\t' + str(dAUC_Perfect) + '\t' + str(dAUC_av_rand)+ '\t' + str(dAUC_D14)+ '\t' +
                str(dAUC_Scorer2)+ '\t' + str(dAUC_Hart17)+ '\t' + str(dAUC_D16woAA)+ '\t' + str(dAUC_D16wAA) + '\t'  +
                str(dAUC_tracrv2_sub) + '\t' + str(dAUC_inDelphi) + '\t' + str(dAUC_Bio_sub) + '\t' + str(dAUC_VBC_sub) +
                '\t' + str(dAUC_tracrv2_rev) + '\t' + str(dAUC_combo_sgRNA) + '\t' + str(dAUC_inDelphi) + '\t' + str(dAUC_Bio_rev) + '\t' + str(dAUC_VBC_rev) +
                '\t' + str(dAUC_CRISPRater) + '\t' + str(dAUC_CRISPRO) + '\t' + str(dAUC_deepCRISPR_nt) + '\t' +
                str(dAUC_deepCRISPR_K562) + '\t' +str(dAUC) + '\t' + str(CEG_depl))
    else:
        line = (data_name + '\t' + str(dAUC_Perfect) + '\t' + str(dAUC_av_rand)+ '\t' + str(dAUC_D14) + '\t' +
                str(dAUC_Scorer2) + '\t' + str(dAUC_Hart17) + '\t' + str(dAUC_D16woAA)+ '\t' + str(dAUC_D16wAA) + '\t' +
                str(dAUC_tracrv2_sub) + '\t' + str(dAUC_inDelphi) + '\t' + str(dAUC_Bio_sub) + '\t' + str(dAUC_VBC_sub) +
                '\t' + str(dAUC_tracrv2_rev) + '\t' + str(dAUC_combo_sgRNA) + '\t' + str(dAUC_inDelphi) + '\t' + str(dAUC_Bio_rev) + '\t' + str(dAUC_VBC_rev) +
                '\t' + str(dAUC_CRISPRater) + '\t' + str('n.a.') + '\t' + str(dAUC_deepCRISPR_nt) + '\t' +
                str('n.a.') + '\t' + str(dAUC) + '\t' + str(CEG_depl))

    outfile_dAUC_best_worst.write('\n' + line)

    ############################# For EVERY SCORE MOUNTED.....

    for number_score,score in enumerate(score_list):
        score_name = score_list_names[number_score]

    ############################# CALCULATE PEARSON CORRELATION MATIRX (every screen with every score)

        Pears, pval = Pearson_plot_score_dependency_lfc(screen_d=Screen_data_d, score_d=score, make_plot=False,
                                          title='lfc'+data_name+'_'+score_name,savepltfolder=plot_folder)
        if number_score == 0:
            pearson_r_lfc_file.write('\n' + data_name + '\t' + str(Pears))
        else:
            pearson_r_lfc_file.write('\t' + str(Pears))
        Pearson_plot_score_dependency_linear(screen_d=Screen_data_d, score_d=score, make_plot=False,
                                          title='lin'+data_name+'_'+score_name,savepltfolder=plot_folder)
        if number_score == 0:
            pearson_r_lin_file.write('\n' + data_name + '\t' + str(Pears))
        else:
            pearson_r_lin_file.write('\t' + str(Pears))

    ############################# CALCULATE dAUC (ROC) (for core essential genes use only best or worst half
    #############################             of sgRNAs as determined by each different scores

        CEG_score, CEG_Rand, CEG_Perfect, nEG_VBC, nEG_Rand, nEG_Perfect = split_screen_by_score(Screen_data_all_d, score)
        dAUC_best,AUC_CEG,AUC_nonEG,sum_total,sum_CEG,sum_nonEG,x,y_dia,y_CEG,y_nonEG = det_dAUC(Screen_data_all_d,CEG_score['best'],nonEG_d)
        dAUC_worst,AUC_CEG,AUC_nonEG,sum_total,sum_CEG,sum_nonEG,x,y_dia,y_CEG,y_nonEG = det_dAUC(Screen_data_all_d,CEG_score['worst'],nonEG_d)
        perfect_dAUC_best,AUC_CEG,AUC_nonEG,sum_total,sum_CEG,sum_nonEG,x,y_dia,y_CEG,y_nonEG = det_dAUC(Screen_data_all_d,CEG_Perfect['best'],nonEG_d)
        perfect_dAUC_worst,AUC_CEG,AUC_nonEG,sum_total,sum_CEG,sum_nonEG,x,y_dia,y_CEG,y_nonEG = det_dAUC(Screen_data_all_d,CEG_Perfect['worst'],nonEG_d)

        ddAUC_perfect = perfect_dAUC_best-perfect_dAUC_worst
        ddAUC = dAUC_best-dAUC_worst
        ddAUC_vs_perfect = ddAUC/ddAUC_perfect
        if number_score == 0:
            ddAUC_output.write('\n' + data_name + '\t' + str(ddAUC))
            rel_to_perfect_ddAUC_output.write('\n' + data_name + '\t' + str(ddAUC_vs_perfect))
        else:
            ddAUC_output.write('\t' + str(ddAUC))
            rel_to_perfect_ddAUC_output.write('\t' + str(ddAUC_vs_perfect))

    x.append(CEG_depl)
    y_VBC.append(dAUC_VBC_rev)
    c_VBC.append('c')
    y_Bio.append(dAUC_Bio_rev)
    c_Bio.append('b')
    y_inD.append(dAUC_inDelphi)
    c_inD.append('orange')
    y_D16.append(dAUC_D16wAA)
    c_D16.append('r')
    y_D16wo.append(dAUC_D16woAA)
    c_D16wo.append('m')
    y_TR2.append(dAUC_Hart17)
    c_TR2.append('r')
    y_rand.append(dAUC_av_rand)

outfile_dAUC_best_worst.close()
pearson_r_lfc_file.close()
pearson_r_lin_file.close()
ddAUC_output.close()

'''
plt.scatter(x,y_Bio,)
plt.scatter(x,y_TR2,)
plt.scatter(x,y_inD)
plt.scatter(x,y_rand)
plt.scatter(x,y_D16wo)
plt.savefig(Analysis_name +'all' '.pdf')
plt.close()
plt.close()
plt.close()
plt.close()
plt.close()
plt.scatter((x,y_VBC))
plt.scatter((x,y_D16))
plt.savefig(Analysis_name +'VBC_vs_ruleset2' '.pdf')
plt.close()
plt.close()
plt.scatter((x,y_VBC))
plt.scatter((x,y_TR2))
plt.savefig(Analysis_name +'VBC_vs_TR2' '.pdf')
'''