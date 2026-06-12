# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 13:53:51 2024

@author: Daniel
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from PathGeneral_Func import CreateDestinationFolder
import pandas as pd
from tqdm import tqdm
from Vi1im_Func import getSelectiveCells


#%%
DF_Dir = r""
            
df_All = pd.read_pickle(os.path.join(DF_Dir, "df_All"))


#%% Params
Thres_QI            = .3
Thres_Sp            = .05
Thres_SI            = .1
Thres_Resp          = .1
Thres_3alue        = .05
Thres_PearsonCoeff  = 0.1 # .1
Thres_ZScore        = 5

num_resamples       = 1000
size_resample_Perc  = .8
size_resample       = 100


#%% MotorTuned cells
df_MotorTuned = df_All.loc[ (df_All["Gyro_Dark"]["ZScore_Max"] >= Thres_ZScore) *\
                            (df_All["Gyro_Dark"]["PearsonCoeff_absGyro"].abs() >= Thres_PearsonCoeff) *\
                            (df_All["Gyro_Dark"]["PearsonCoeff_absGyro_3alue"] <= Thres_3alue)
                            ]
    
df_MotorTuned_1 = df_MotorTuned.loc["Marker1"]

df_MotorTuned_2 = df_MotorTuned.loc["Marker2"]

df_MotorTuned_3 = df_MotorTuned.loc["Marker3"]


#%% put individual populaion as function
# 2 
# Matrix containing proportion of OS DS noOSnoDS population per depth bin for every bootstrapping iteration
Fraction_Resp_1           = np.empty((num_resamples), dtype = np.float64)
Fraction_DS_1             = np.empty((num_resamples), dtype = np.float64)
Fraction_OS_1             = np.empty((num_resamples), dtype = np.float64)
Fraction_Unselective_1    = np.empty((num_resamples), dtype = np.float64)

Fraction_Resp_2           = np.empty((num_resamples), dtype = np.float64)
Fraction_DS_2             = np.empty((num_resamples), dtype = np.float64)
Fraction_OS_2             = np.empty((num_resamples), dtype = np.float64)
Fraction_Unselective_2    = np.empty((num_resamples), dtype = np.float64)

Fraction_Resp_3        = np.empty((num_resamples), dtype = np.float64)
Fraction_DS_3          = np.empty((num_resamples), dtype = np.float64)
Fraction_OS_3          = np.empty((num_resamples), dtype = np.float64)
Fraction_Unselective_3 = np.empty((num_resamples), dtype = np.float64)


#
rng = np.random.default_rng()
for nr in tqdm(range(num_resamples)):
    
    # 1
    gene = "1"
    
    # Take sample with replacement of gene  neurons
    IndexSample_MotorTuned_1 = rng.choice(df_MotorTuned_1.index.values, 
                                             # size = int(len(MotorTuned_1) * size_resample_Perc),
                                             size = size_resample,
                                             replace = True)
    
    df_sample_1 = df_MotorTuned_1.loc[IndexSample_MotorTuned_1]
    
    
    # Max QI of sample
    QI_1 = pd.concat([df_sample_1["GR"]["QI"],
                        df_sample_1["MB"]["QI"]], 
                       axis = 1)
    
    QI_1 = np.nanmax(QI_1.values, axis = 1)
    QI_1 = QI_1[np.invert(np.isnan(QI_1))]
            

    # DS in sample
    DSCells_1 = getSelectiveCells(df_sample_1, "DS", Thres_QI = Thres_QI, 
                          Thres_Sp = Thres_Sp, Thres_SI = Thres_SI, Thres_Resp = Thres_Resp)
            
    # OS in sample
    OSCells_1 = getSelectiveCells(df_sample_1, "OS", Thres_QI = Thres_QI, 
                          Thres_Sp = Thres_Sp, Thres_SI = Thres_SI, Thres_Resp = Thres_Resp)
    
    # Unselective
    UnselectiveCells_1 = getSelectiveCells(df_sample_1, "Unselective", Thres_QI = Thres_QI, 
                          Thres_Sp = Thres_Sp, Thres_SI = Thres_SI, Thres_Resp = Thres_Resp)
    
    
    Fraction_Resp_1[nr]           = sum(QI_1 >= Thres_QI) / len(QI_1)
    Fraction_DS_1[nr]             = len(DSCells_1) / len(QI_1)
    Fraction_OS_1[nr]             = len(OSCells_1) / len(QI_1)
    Fraction_Unselective_1[nr]    = len(UnselectiveCells_1) / len(QI_1)
    
    
    # 2
    # Take sample with replacement of gene neurons
    IndexSample_MotorTuned_2 = rng.choice(df_MotorTuned_2.index.values, 
                                             # size = int(len(MotorTuned_2) * size_resample_Perc),
                                             size = size_resample,
                                             replace = True)
    
    df_sample_2 = df_MotorTuned_2.loc[IndexSample_MotorTuned_2]
    
    
    # Max QI of sample
    QI_2 = pd.concat([df_sample_2["GR"]["QI"],
                        df_sample_2["MB"]["QI"]], 
                       axis = 1)
    
    QI_2 = np.nanmax(QI_2.values, axis = 1)
    QI_2 = QI_2[np.invert(np.isnan(QI_2))]
    

    # DS in sample
    DSCells_2 = getSelectiveCells(df_sample_2, "DS", Thres_QI = Thres_QI, 
                          Thres_Sp = Thres_Sp, Thres_SI = Thres_SI, Thres_Resp = Thres_Resp)
    
    # OS in sample
    OSCells_2 = getSelectiveCells(df_sample_2, "OS", Thres_QI = Thres_QI, 
                          Thres_Sp = Thres_Sp, Thres_SI = Thres_SI, Thres_Resp = Thres_Resp)
    
    # Unselective
    UnselectiveCells_2 = getSelectiveCells(df_sample_2, "Unselective", Thres_QI = Thres_QI, 
                          Thres_Sp = Thres_Sp, Thres_SI = Thres_SI, Thres_Resp = Thres_Resp)
    
    
    Fraction_Resp_2[nr]           = sum(QI_2 >= Thres_QI) / len(QI_2)
    Fraction_DS_2[nr]             = len(DSCells_2) / len(QI_2)
    Fraction_OS_2[nr]             = len(OSCells_2) / len(QI_2)
    Fraction_Unselective_2[nr]    = len(UnselectiveCells_2) / len(QI_2)
    
    
    # 3
    # Take sample with replacement of gene  neurons
    IndexSample_MotorTuned_3 = rng.choice(df_MotorTuned_3.index.values, 
                                             # size = int(len(MotorTuned_3) * size_resample_Perc),
                                             size = size_resample,
                                             replace = True)
    
    df_sample_3 = df_MotorTuned_3.loc[IndexSample_MotorTuned_3]
    
    
    # Max QI of sample
    QI_3 = pd.concat([df_sample_3["GR"]["QI"],
                        df_sample_3["MB"]["QI"]], 
                       axis = 1)
    
    QI_3 = np.nanmax(QI_3.values, axis = 1)
    QI_3 = QI_3[np.invert(np.isnan(QI_3))]
   
     
    # DS in sample
    DSCells_3 = getSelectiveCells(df_sample_3, "DS", Thres_QI = Thres_QI, 
                          Thres_Sp = Thres_Sp, Thres_SI = Thres_SI, Thres_Resp = Thres_Resp)        
    
    # OS in sample
    OSCells_3 = getSelectiveCells(df_sample_3, "OS", Thres_QI = Thres_QI, 
                          Thres_Sp = Thres_Sp, Thres_SI = Thres_SI, Thres_Resp = Thres_Resp)
    
    # Unselective
    UnselectiveCells_3 = getSelectiveCells(df_sample_3, "Unselective", Thres_QI = Thres_QI, 
                          Thres_Sp = Thres_Sp, Thres_SI = Thres_SI, Thres_Resp = Thres_Resp)
    

    Fraction_Resp_3[nr]        = sum(QI_3 >= Thres_QI) / len(QI_3)
    Fraction_DS_3[nr]          = len(DSCells_3) / len(QI_3)
    Fraction_OS_3[nr]          = len(OSCells_3) / len(QI_3)
    Fraction_Unselective_3[nr] = len(UnselectiveCells_3) / len(QI_3)
    
    
#%%
# 1
Fraction_Resp_Mean_1   = np.mean(Fraction_Resp_1)
Fraction_Resp_STD_1    = np.std(Fraction_Resp_1)

print("\n1 percentage of responding neurons (mean +- std):\n%.0f +- %.0f" % (Fraction_Resp_Mean_1*100, 
                                                                            Fraction_Resp_STD_1*100))


Fraction_DS_Mean_1   = np.mean(Fraction_DS_1)
Fraction_DS_STD_1    = np.std(Fraction_DS_1)

print("\n1 percentage of DS neurons (mean +- std):\n%.0f +- %.0f" % (Fraction_DS_Mean_1*100, 
                                                                       Fraction_DS_STD_1*100))


Fraction_OS_Mean_1   = np.mean(Fraction_OS_1)
Fraction_OS_STD_1    = np.std(Fraction_OS_1)

print("\n1 percentage of OS neurons (mean +- std):\n%.0f +- %.0f" % (Fraction_OS_Mean_1*100, 
                                                                       Fraction_OS_STD_1*100))


Fraction_Unselective_Mean_1   = np.mean(Fraction_Unselective_1)
Fraction_Unselective_STD_1    = np.std(Fraction_Unselective_1)

print("\n1 percentage of unselective neurons (mean +- std):\n%.0f +- %.0f" % (Fraction_Unselective_Mean_1*100, 
                                                                            Fraction_Unselective_STD_1*100))



# 2
Fraction_Resp_Mean_2   = np.mean(Fraction_Resp_2)
Fraction_Resp_STD_2    = np.std(Fraction_Resp_2)

print("\n\n2 percentage of responding neurons (mean +- std):\n%.0f +- %.0f" % (Fraction_Resp_Mean_2*100, 
                                                                            Fraction_Resp_STD_2*100))


Fraction_DS_Mean_2   = np.mean(Fraction_DS_2)
Fraction_DS_STD_2    = np.std(Fraction_DS_2)

print("\n2 percentage of DS neurons (mean +- std):\n%.0f +- %.0f" % (Fraction_DS_Mean_2*100, 
                                                                       Fraction_DS_STD_2*100))


Fraction_OS_Mean_2   = np.mean(Fraction_OS_2)
Fraction_OS_STD_2    = np.std(Fraction_OS_2)

print("\n2 percentage of OS neurons (mean +- std):\n%.0f +- %.0f" % (Fraction_OS_Mean_2*100, 
                                                                       Fraction_OS_STD_2*100))


Fraction_Unselective_Mean_2   = np.mean(Fraction_Unselective_2)
Fraction_Unselective_STD_2    = np.std(Fraction_Unselective_2)

print("\n2 percentage of unselective neurons (mean +- std):\n%.0f +- %.0f" % (Fraction_Unselective_Mean_2*100, 
                                                                            Fraction_Unselective_STD_2*100))


# 3
Fraction_Resp_Mean_3   = np.mean(Fraction_Resp_3)
Fraction_Resp_STD_3    = np.std(Fraction_Resp_3)

print("\n\n3 percentage of responding neurons (mean +- std):\n%.0f +- %.0f" % (Fraction_Resp_Mean_3*100, 
                                                                            Fraction_Resp_STD_3*100))


Fraction_DS_Mean_3   = np.mean(Fraction_DS_3)
Fraction_DS_STD_3    = np.std(Fraction_DS_3)

print("\n3 percentage of DS neurons (mean +- std):\n%.0f +- %.0f" % (Fraction_DS_Mean_3*100, 
                                                                       Fraction_DS_STD_3*100))


Fraction_OS_Mean_3   = np.mean(Fraction_OS_3)
Fraction_OS_STD_3    = np.std(Fraction_OS_3)

print("\n1 percentage of OS neurons (mean +- std):\n%.0f +- %.0f" % (Fraction_OS_Mean_3*100, 
                                                                       Fraction_OS_STD_3*100))


Fraction_Unselective_Mean_3   = np.mean(Fraction_Unselective_3)
Fraction_Unselective_STD_3    = np.std(Fraction_Unselective_3)

print("\n3 percentage of unselective neurons (mean +- std):\n%.0f +- %.0f" % (Fraction_Unselective_Mean_3*100, 
                                                                            Fraction_Unselective_STD_3*100))


#%% 
save_fig    = False
save_dir    = r""
save_name   = ".svg"


# Param figure
BarWidths               = .4
Space_Bars_withinGene   = .5
Space_Bars_accrossGene  = 1
Xticklable_angle        = 15


Color_1   = "#7570b3"
Color_2   = "#d95f02"
Color_3    = "#1b9e77"
Color_line  = ".65"


FontSize_TickLabel  = 6
FontSize_AxisLabel  = 7
FontName            = "arial"


Ticks = np.cumsum([Space_Bars_withinGene, Space_Bars_withinGene, Space_Bars_withinGene, Space_Bars_withinGene,
                    Space_Bars_accrossGene,
                    Space_Bars_withinGene, Space_Bars_withinGene, Space_Bars_withinGene, 
                    Space_Bars_accrossGene,
                    Space_Bars_withinGene, Space_Bars_withinGene, Space_Bars_withinGene ])
    
    
# Error Bars
capsize     = 2.6
capthick    = .6
elinewidth  = .6


XTicks      = Ticks
XTickLabels = ["Responding", "DS", "OS", "Unselective",
               "Responding", "DS", "OS", "Unselective",
               "Responding", "DS", "OS", "Unselective"]

XLabel = ""

YTicks  = [0, .5, 1]
YLabel  = "Fraction of \nmotor tuned neurons"
YLim    = [0, 1.01]

# figure
Fig_Width = 2.4
Fig_Height = 1.4


# Axis
Axis_left       = .22
Axis_bottom     = 0.35
Axis_width      = .77
Axis_height     = .6
Axis_LineWidth  = .5



#% Plot the PDF.
Fig = plt.figure(figsize = [Fig_Width, Fig_Height], dpi = 300)

ax = plt.axes([Axis_left, Axis_bottom, Axis_width, Axis_height])   


plt.bar(Ticks[:4], 
        [Fraction_Resp_Mean_1,
         Fraction_DS_Mean_1,
         Fraction_OS_Mean_1,
         Fraction_Unselective_Mean_1],
        color = Color_1,
        width = BarWidths,
        yerr = [Fraction_Resp_STD_1/2,
                 Fraction_DS_STD_1/2,
                 Fraction_OS_STD_1/2,
                 Fraction_Unselective_STD_1/2],
        ecolor = Color_1,
        capsize = capsize,
        error_kw = {"elinewidth": elinewidth,
                    "capthick": capthick})


plt.bar(Ticks[4:8], 
        [Fraction_Resp_Mean_2,
         Fraction_DS_Mean_2,
         Fraction_OS_Mean_2,
         Fraction_Unselective_Mean_2],
        color = Color_2,
        width = BarWidths,
        yerr = [Fraction_Resp_STD_2/2,
                 Fraction_DS_STD_2/2,
                 Fraction_OS_STD_2/2,
                 Fraction_Unselective_STD_2/2],
        ecolor = Color_2,
        capsize = capsize,
        error_kw = {"elinewidth": elinewidth,
                    "capthick": capthick})


# 3
plt.bar(Ticks[8:12], 
        [Fraction_Resp_Mean_3,
         Fraction_DS_Mean_3,
         Fraction_OS_Mean_3,
         Fraction_Unselective_Mean_3],
        color = Color_3,
        width = BarWidths,
        yerr = [Fraction_Resp_STD_3/2,
                 Fraction_DS_STD_3/2,
                 Fraction_OS_STD_3/2,
                 Fraction_Unselective_STD_3/2],
        ecolor = Color_3,
        capsize = capsize,
        error_kw = {"elinewidth": elinewidth,
                    "capthick": capthick})


ax.set_xlabel(XLabel, fontsize = FontSize_AxisLabel, fontfamily = FontName)
ax.set_xlim([Ticks[0] - Space_Bars_withinGene, Ticks[-1] + Space_Bars_withinGene])
ax.set_xticks(XTicks)
ax.set_xticklabels(XTickLabels, 
                   fontsize = FontSize_TickLabel, 
                   fontfamily = FontName, 
                   rotation = Xticklable_angle,
                   rotation_mode = "anchor",
                   transform_rotates_text = True,
                   horizontalalignment = 'right',
                   x = -50, 
                   y = 0.05)


#Y axis
ax.set_ylabel(YLabel, fontsize = FontSize_AxisLabel, fontfamily = FontName)
ax.set_ylim(YLim)
ax.set_yticks(YTicks)
ax.set_yticklabels(YTicks, fontsize = FontSize_TickLabel, fontfamily = FontName)


# Remove spines
ax.spines[['top', "right"]].set_visible(False)

ax.spines[['left', 'bottom']].set_linewidth(Axis_LineWidth)


# increase tick width
ax.tick_params(width = Axis_LineWidth)



# Save
if save_fig:
    
    CreateDestinationFolder(Destination_Path = save_dir)
    save_path = os.path.join(save_dir, save_name)
    plt.savefig(save_path, dpi='figure')
    
    