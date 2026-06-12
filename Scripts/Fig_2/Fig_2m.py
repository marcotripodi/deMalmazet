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
from VisStim_Func import getSelectiveCells


#%%
DF_Dir = r"O:\OneDrive - MRC Laboratory of Molecular Biology\Work\Projects\Ongoing\Inhibition\Experiments\Completed\Imaging\All\Code\MB_GR"
            
df_All = pd.read_pickle(os.path.join(DF_Dir, "df_All"))


#%% Params
Thres_QI            = .3
Thres_Sp            = .05
Thres_SI            = .1
Thres_Resp          = .1
Thres_PValue        = .05
Thres_PearsonCoeff  = 0.1 # .1
Thres_ZScore        = 5

num_resamples       = 1000
size_resample_Perc  = .8
size_resample       = 100


#%% MotorTuned cells
df_MotorTuned = df_All.loc[ (df_All["Gyro_Dark"]["ZScore_Max"] >= Thres_ZScore) *\
                            (df_All["Gyro_Dark"]["PearsonCoeff_absGyro"].abs() >= Thres_PearsonCoeff) *\
                            (df_All["Gyro_Dark"]["PearsonCoeff_absGyro_PValue"] <= Thres_PValue)
                            ]
    
df_MotorTuned_SST = df_MotorTuned.loc["SST"]

df_MotorTuned_CCK = df_MotorTuned.loc["CCK"]

df_MotorTuned_PV = df_MotorTuned.loc["PV"]


#%% put individual populaion as function
# 2 
# Matrix containing proportion of OS DS noOSnoDS population per depth bin for every bootstrapping iteration
Fraction_Resp_SST           = np.empty((num_resamples), dtype = np.float64)
Fraction_DS_SST             = np.empty((num_resamples), dtype = np.float64)
Fraction_OS_SST             = np.empty((num_resamples), dtype = np.float64)
Fraction_Unselective_SST    = np.empty((num_resamples), dtype = np.float64)

Fraction_Resp_CCK           = np.empty((num_resamples), dtype = np.float64)
Fraction_DS_CCK             = np.empty((num_resamples), dtype = np.float64)
Fraction_OS_CCK             = np.empty((num_resamples), dtype = np.float64)
Fraction_Unselective_CCK    = np.empty((num_resamples), dtype = np.float64)

Fraction_Resp_PV        = np.empty((num_resamples), dtype = np.float64)
Fraction_DS_PV          = np.empty((num_resamples), dtype = np.float64)
Fraction_OS_PV          = np.empty((num_resamples), dtype = np.float64)
Fraction_Unselective_PV = np.empty((num_resamples), dtype = np.float64)


#
rng = np.random.default_rng()
for nr in tqdm(range(num_resamples)):
    
    # SST
    gene = "SST"
    
    # Take sample with replacement of gene  neurons
    IndexSample_MotorTuned_SST = rng.choice(df_MotorTuned_SST.index.values, 
                                             # size = int(len(MotorTuned_SST) * size_resample_Perc),
                                             size = size_resample,
                                             replace = True)
    
    df_sample_SST = df_MotorTuned_SST.loc[IndexSample_MotorTuned_SST]
    
    
    # Max QI of sample
    QI_SST = pd.concat([df_sample_SST["GR"]["QI"],
                        df_sample_SST["MB"]["QI"]], 
                       axis = 1)
    
    QI_SST = np.nanmax(QI_SST.values, axis = 1)
    QI_SST = QI_SST[np.invert(np.isnan(QI_SST))]
            

    # DS in sample
    DSCells_SST = getSelectiveCells(df_sample_SST, "DS", Thres_QI = Thres_QI, 
                          Thres_Sp = Thres_Sp, Thres_SI = Thres_SI, Thres_Resp = Thres_Resp)
            
    # OS in sample
    OSCells_SST = getSelectiveCells(df_sample_SST, "OS", Thres_QI = Thres_QI, 
                          Thres_Sp = Thres_Sp, Thres_SI = Thres_SI, Thres_Resp = Thres_Resp)
    
    # Unselective
    UnselectiveCells_SST = getSelectiveCells(df_sample_SST, "Unselective", Thres_QI = Thres_QI, 
                          Thres_Sp = Thres_Sp, Thres_SI = Thres_SI, Thres_Resp = Thres_Resp)
    
    
    Fraction_Resp_SST[nr]           = sum(QI_SST >= Thres_QI) / len(QI_SST)
    Fraction_DS_SST[nr]             = len(DSCells_SST) / len(QI_SST)
    Fraction_OS_SST[nr]             = len(OSCells_SST) / len(QI_SST)
    Fraction_Unselective_SST[nr]    = len(UnselectiveCells_SST) / len(QI_SST)
    
    
    # CCK
    # Take sample with replacement of gene neurons
    IndexSample_MotorTuned_CCK = rng.choice(df_MotorTuned_CCK.index.values, 
                                             # size = int(len(MotorTuned_CCK) * size_resample_Perc),
                                             size = size_resample,
                                             replace = True)
    
    df_sample_CCK = df_MotorTuned_CCK.loc[IndexSample_MotorTuned_CCK]
    
    
    # Max QI of sample
    QI_CCK = pd.concat([df_sample_CCK["GR"]["QI"],
                        df_sample_CCK["MB"]["QI"]], 
                       axis = 1)
    
    QI_CCK = np.nanmax(QI_CCK.values, axis = 1)
    QI_CCK = QI_CCK[np.invert(np.isnan(QI_CCK))]
    

    # DS in sample
    DSCells_CCK = getSelectiveCells(df_sample_CCK, "DS", Thres_QI = Thres_QI, 
                          Thres_Sp = Thres_Sp, Thres_SI = Thres_SI, Thres_Resp = Thres_Resp)
    
    # OS in sample
    OSCells_CCK = getSelectiveCells(df_sample_CCK, "OS", Thres_QI = Thres_QI, 
                          Thres_Sp = Thres_Sp, Thres_SI = Thres_SI, Thres_Resp = Thres_Resp)
    
    # Unselective
    UnselectiveCells_CCK = getSelectiveCells(df_sample_CCK, "Unselective", Thres_QI = Thres_QI, 
                          Thres_Sp = Thres_Sp, Thres_SI = Thres_SI, Thres_Resp = Thres_Resp)
    
    
    Fraction_Resp_CCK[nr]           = sum(QI_CCK >= Thres_QI) / len(QI_CCK)
    Fraction_DS_CCK[nr]             = len(DSCells_CCK) / len(QI_CCK)
    Fraction_OS_CCK[nr]             = len(OSCells_CCK) / len(QI_CCK)
    Fraction_Unselective_CCK[nr]    = len(UnselectiveCells_CCK) / len(QI_CCK)
    
    
    # PV
    # Take sample with replacement of gene  neurons
    IndexSample_MotorTuned_PV = rng.choice(df_MotorTuned_PV.index.values, 
                                             # size = int(len(MotorTuned_PV) * size_resample_Perc),
                                             size = size_resample,
                                             replace = True)
    
    df_sample_PV = df_MotorTuned_PV.loc[IndexSample_MotorTuned_PV]
    
    
    # Max QI of sample
    QI_PV = pd.concat([df_sample_PV["GR"]["QI"],
                        df_sample_PV["MB"]["QI"]], 
                       axis = 1)
    
    QI_PV = np.nanmax(QI_PV.values, axis = 1)
    QI_PV = QI_PV[np.invert(np.isnan(QI_PV))]
   
     
    # DS in sample
    DSCells_PV = getSelectiveCells(df_sample_PV, "DS", Thres_QI = Thres_QI, 
                          Thres_Sp = Thres_Sp, Thres_SI = Thres_SI, Thres_Resp = Thres_Resp)        
    
    # OS in sample
    OSCells_PV = getSelectiveCells(df_sample_PV, "OS", Thres_QI = Thres_QI, 
                          Thres_Sp = Thres_Sp, Thres_SI = Thres_SI, Thres_Resp = Thres_Resp)
    
    # Unselective
    UnselectiveCells_PV = getSelectiveCells(df_sample_PV, "Unselective", Thres_QI = Thres_QI, 
                          Thres_Sp = Thres_Sp, Thres_SI = Thres_SI, Thres_Resp = Thres_Resp)
    

    Fraction_Resp_PV[nr]        = sum(QI_PV >= Thres_QI) / len(QI_PV)
    Fraction_DS_PV[nr]          = len(DSCells_PV) / len(QI_PV)
    Fraction_OS_PV[nr]          = len(OSCells_PV) / len(QI_PV)
    Fraction_Unselective_PV[nr] = len(UnselectiveCells_PV) / len(QI_PV)
    
    
#%%
# SST
Fraction_Resp_Mean_SST   = np.mean(Fraction_Resp_SST)
Fraction_Resp_STD_SST    = np.std(Fraction_Resp_SST)

print("\nSST percentage of responding neurons (mean +- std):\n%.0f +- %.0f" % (Fraction_Resp_Mean_SST*100, 
                                                                            Fraction_Resp_STD_SST*100))


Fraction_DS_Mean_SST   = np.mean(Fraction_DS_SST)
Fraction_DS_STD_SST    = np.std(Fraction_DS_SST)

print("\nSST percentage of DS neurons (mean +- std):\n%.0f +- %.0f" % (Fraction_DS_Mean_SST*100, 
                                                                       Fraction_DS_STD_SST*100))


Fraction_OS_Mean_SST   = np.mean(Fraction_OS_SST)
Fraction_OS_STD_SST    = np.std(Fraction_OS_SST)

print("\nSST percentage of OS neurons (mean +- std):\n%.0f +- %.0f" % (Fraction_OS_Mean_SST*100, 
                                                                       Fraction_OS_STD_SST*100))


Fraction_Unselective_Mean_SST   = np.mean(Fraction_Unselective_SST)
Fraction_Unselective_STD_SST    = np.std(Fraction_Unselective_SST)

print("\nSST percentage of unselective neurons (mean +- std):\n%.0f +- %.0f" % (Fraction_Unselective_Mean_SST*100, 
                                                                            Fraction_Unselective_STD_SST*100))



# CCK
Fraction_Resp_Mean_CCK   = np.mean(Fraction_Resp_CCK)
Fraction_Resp_STD_CCK    = np.std(Fraction_Resp_CCK)

print("\n\nCCK percentage of responding neurons (mean +- std):\n%.0f +- %.0f" % (Fraction_Resp_Mean_CCK*100, 
                                                                            Fraction_Resp_STD_CCK*100))


Fraction_DS_Mean_CCK   = np.mean(Fraction_DS_CCK)
Fraction_DS_STD_CCK    = np.std(Fraction_DS_CCK)

print("\nCCK percentage of DS neurons (mean +- std):\n%.0f +- %.0f" % (Fraction_DS_Mean_CCK*100, 
                                                                       Fraction_DS_STD_CCK*100))


Fraction_OS_Mean_CCK   = np.mean(Fraction_OS_CCK)
Fraction_OS_STD_CCK    = np.std(Fraction_OS_CCK)

print("\nCCK percentage of OS neurons (mean +- std):\n%.0f +- %.0f" % (Fraction_OS_Mean_CCK*100, 
                                                                       Fraction_OS_STD_CCK*100))


Fraction_Unselective_Mean_CCK   = np.mean(Fraction_Unselective_CCK)
Fraction_Unselective_STD_CCK    = np.std(Fraction_Unselective_CCK)

print("\nCCK percentage of unselective neurons (mean +- std):\n%.0f +- %.0f" % (Fraction_Unselective_Mean_CCK*100, 
                                                                            Fraction_Unselective_STD_CCK*100))


# PV
Fraction_Resp_Mean_PV   = np.mean(Fraction_Resp_PV)
Fraction_Resp_STD_PV    = np.std(Fraction_Resp_PV)

print("\n\nPV percentage of responding neurons (mean +- std):\n%.0f +- %.0f" % (Fraction_Resp_Mean_PV*100, 
                                                                            Fraction_Resp_STD_PV*100))


Fraction_DS_Mean_PV   = np.mean(Fraction_DS_PV)
Fraction_DS_STD_PV    = np.std(Fraction_DS_PV)

print("\nPV percentage of DS neurons (mean +- std):\n%.0f +- %.0f" % (Fraction_DS_Mean_PV*100, 
                                                                       Fraction_DS_STD_PV*100))


Fraction_OS_Mean_PV   = np.mean(Fraction_OS_PV)
Fraction_OS_STD_PV    = np.std(Fraction_OS_PV)

print("\nSST percentage of OS neurons (mean +- std):\n%.0f +- %.0f" % (Fraction_OS_Mean_PV*100, 
                                                                       Fraction_OS_STD_PV*100))


Fraction_Unselective_Mean_PV   = np.mean(Fraction_Unselective_PV)
Fraction_Unselective_STD_PV    = np.std(Fraction_Unselective_PV)

print("\nPV percentage of unselective neurons (mean +- std):\n%.0f +- %.0f" % (Fraction_Unselective_Mean_PV*100, 
                                                                            Fraction_Unselective_STD_PV*100))


#%% 
save_fig    = False
save_dir    = r"O:\OneDrive - MRC Laboratory of Molecular Biology\General - VGAT Paper\Figures\Paper\Figure4\Ver5\BarPlot"
save_name   = "BarPlot_MotorTuned_WithErrorBars_Resp_unselective_Selective.svg"


# Param figure
BarWidths               = .4
Space_Bars_withinGene   = .5
Space_Bars_accrossGene  = 1
Xticklable_angle        = 15


Color_SST   = "#7570b3"
Color_CCK   = "#d95f02"
Color_PV    = "#1b9e77"
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
        [Fraction_Resp_Mean_SST,
         Fraction_DS_Mean_SST,
         Fraction_OS_Mean_SST,
         Fraction_Unselective_Mean_SST],
        color = Color_SST,
        width = BarWidths,
        yerr = [Fraction_Resp_STD_SST/2,
                 Fraction_DS_STD_SST/2,
                 Fraction_OS_STD_SST/2,
                 Fraction_Unselective_STD_SST/2],
        ecolor = Color_SST,
        capsize = capsize,
        error_kw = {"elinewidth": elinewidth,
                    "capthick": capthick})


plt.bar(Ticks[4:8], 
        [Fraction_Resp_Mean_CCK,
         Fraction_DS_Mean_CCK,
         Fraction_OS_Mean_CCK,
         Fraction_Unselective_Mean_CCK],
        color = Color_CCK,
        width = BarWidths,
        yerr = [Fraction_Resp_STD_CCK/2,
                 Fraction_DS_STD_CCK/2,
                 Fraction_OS_STD_CCK/2,
                 Fraction_Unselective_STD_CCK/2],
        ecolor = Color_CCK,
        capsize = capsize,
        error_kw = {"elinewidth": elinewidth,
                    "capthick": capthick})


# PV
plt.bar(Ticks[8:12], 
        [Fraction_Resp_Mean_PV,
         Fraction_DS_Mean_PV,
         Fraction_OS_Mean_PV,
         Fraction_Unselective_Mean_PV],
        color = Color_PV,
        width = BarWidths,
        yerr = [Fraction_Resp_STD_PV/2,
                 Fraction_DS_STD_PV/2,
                 Fraction_OS_STD_PV/2,
                 Fraction_Unselective_STD_PV/2],
        ecolor = Color_PV,
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
    
    