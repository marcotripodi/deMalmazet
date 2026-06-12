# -*- coding: utf-8 -*-
"""
Created on Thu Jun 19 14:57:14 2025

@author: Daniel
"""

import pickle
import os
import matplotlib.pyplot as plt
import numpy as np
from PathGeneral_Func import CreateDestinationFolder

from scipy import stats
from Hits_Func import GetHistContrast_SingleMouse


#%% Load Data
# All mice
Hits_Save_Dir = r"W:\Private_Lab\Daniel\MovingSpotArena\VGAT-Cre x Sst-Flpo\SingleDot_LargeExpression"


Hits_Save_Name = "HitTimes_SST_SingleDot"


# Load hit counts
with open(os.path.join(Hits_Save_Dir, Hits_Save_Name), 'rb') as handle:
    Hits_Time_Contrasts_AllMice = pickle.load(handle)
    
    
Days_DCZ = ['Day3']
Days_PBS = ["Day1", "Day2"]


# Days_PBS1 = ["Day1"]
# Days_PBS2 = ['Day2']


# Days_DCZ = ['Day3']
# Days_PBS = ["Day2"]


Days_PBS_int = [int(day[day.find("y") + 1:]) for day in Days_PBS]
Days_DCZ_int = [int(day[day.find("y") + 1:]) for day in Days_DCZ]

# NStim per contrast
nTrials = 3
nAngles = 8
nStimPerContrast = nTrials * nAngles  

Contrast    = [5, 20, 100]
nContrasts  = len(Contrast)


#%% Get mean Hits per contrast per mouse
nDays_PBS = len(Days_PBS)
nDays_DCZ = len(Days_DCZ)

nMice = len(Hits_Time_Contrasts_AllMice)


Hits_Contrast_All_DCZ = np.empty((nMice, nDays_DCZ))
Hits_Contrast_100_DCZ = np.empty((nMice, nDays_DCZ))
Hits_Contrast_020_DCZ = np.empty((nMice, nDays_DCZ))
Hits_Contrast_005_DCZ = np.empty((nMice, nDays_DCZ))

Hits_Contrast_All_PBS = np.empty((nMice, nDays_PBS))
Hits_Contrast_100_PBS = np.empty((nMice, nDays_PBS))
Hits_Contrast_020_PBS = np.empty((nMice, nDays_PBS))
Hits_Contrast_005_PBS = np.empty((nMice, nDays_PBS))

mouse_num = 0
for mouse in Hits_Time_Contrasts_AllMice:  

    HitCount_DCZ = GetHistContrast_SingleMouse(Hits_Time_Contrasts_AllMice[mouse], Days_DCZ)
    
    
    # Convert to rate
    HitCount_DCZ /= nStimPerContrast
    
    
    Hits_Contrast_All_DCZ[mouse_num, :] = HitCount_DCZ.sum(axis = 1) / nContrasts
    Hits_Contrast_100_DCZ[mouse_num, :] = HitCount_DCZ.loc[:, -1.0]
    Hits_Contrast_020_DCZ[mouse_num, :] = HitCount_DCZ.loc[:, HitCount_DCZ.columns[1]]
    Hits_Contrast_005_DCZ[mouse_num, :] = HitCount_DCZ.loc[:, HitCount_DCZ.columns[0]]
    
    
    # PBS
    HitCount_PBS = GetHistContrast_SingleMouse(Hits_Time_Contrasts_AllMice[mouse], Days_PBS)
    
    HitCount_PBS /= nStimPerContrast
    
    Hits_Contrast_All_PBS[mouse_num, :] = HitCount_PBS.sum(axis = 1) / nContrasts
    Hits_Contrast_100_PBS[mouse_num, :] = HitCount_PBS.loc[:, -1.0]
    Hits_Contrast_020_PBS[mouse_num, :] = HitCount_PBS.loc[:, HitCount_PBS.columns[1]]
    Hits_Contrast_005_PBS[mouse_num, :] = HitCount_PBS.loc[:, HitCount_PBS.columns[0]]
    
    mouse_num += 1
    
    
#%% Get data to use
Thres_Perf = .1 #3 / nStimPerContrast
Thres_Perf = .05 #3 / nStimPerContrast


# Data2Use_PBS = Hits_Contrast_All_PBS
Data2Use_PBS = Hits_Contrast_100_PBS.copy()
Data2Use_PBS = Hits_Contrast_020_PBS.copy()
Data2Use_PBS = Hits_Contrast_005_PBS.copy()

# Data2Use_DCZ = Hits_Contrast_All_DCZ
Data2Use_DCZ = Hits_Contrast_100_DCZ.copy()
Data2Use_DCZ = Hits_Contrast_020_DCZ.copy()
Data2Use_DCZ = Hits_Contrast_005_DCZ.copy()


Data2Use_PBS = np.mean(Data2Use_PBS, axis = 1, keepdims = True)

Points2Keep = np.flatnonzero(Data2Use_PBS > Thres_Perf)

Data2Use_PBS = Data2Use_PBS[Points2Keep, :]
Data2Use_DCZ = Data2Use_DCZ[Points2Keep, :]


#%% All
Data2Use = np.concat((Data2Use_PBS,
                      Data2Use_DCZ),
                     axis = 1)
                      
Days2Use_int = Days_PBS_int #+ Days_DCZ_int


# Sort by days
IdxSort = np.argsort(Days2Use_int)
Data2Use = Data2Use[:, IdxSort]
Days2Use_int = np.sort(Days2Use_int)


#%%

# Save
Fig_SaveFlag = False
Fig_SaveDir = r"O:\OneDrive - MRC Laboratory of Molecular Biology\General - VGAT Paper\Figures\Paper\Figure6\Ver6\Perf"

# Fig_SaveName = "Peformance_over_days_All.svg"
# Fig_SaveName = "Peformance_over_days_100.svg"
Fig_SaveName = "Peformance_over_days_020.svg"
# Fig_SaveName = "Peformance_over_days_005.svg"


# figure
Fig_Width   = 1.2
Fig_Height  = 1.5


# Axis
Axis_left       = .28
Axis_bottom     = 0.25
Axis_width      = .65
Axis_height     = .6
Axis_LineWidth  = .5
    
# Legend
Legend_ShowFlag         = False
Legend_X0               = .1
Legend_Y0               = .8
Legend_ColumnSpacing    = .4
Legend_handlelength     = .8


# Mean acrross mice
MeanErrorBars_LineWidth     = .9
MeanErrorBars_Color = "0.1"


# Single mice lines
SingleMice_LineWidth     = .4
SingleMice_Color        = ".4"
# SingleMice_Color    = None
SingleMice_Marker       = "."  # "x"
SingleMice_MarkerSize   = 2


FontSize_TickLabel  = 6
FontSize_AxisLabel  = 7
FontName            = "arial"


# Y axis
YMax        = np.max(Data2Use)
YTicks      = np.arange(0, YMax, step = .1)
# YTicks      = np.arange(0, YMax, step = .05)
YTickLabels = (YTicks*100).astype(np.int32)
YLabel      = "% dot interceptions"
# YLabel      = "% interceptions per dot presentation"
YSlack      = .01
YLim        = [0, YMax + YSlack]



# Xaxis
XSlack      = .2
X_Lim       = [Days2Use_int[0] - XSlack, Days2Use_int[-1] + XSlack]
XTicks      = np.arange(Days2Use_int[0], Days2Use_int[-1]+1, step = 1) #[1, 5, 10, 15]

# XTickLabels = XTicks - 1
# XLabel      = "Days"

XTickLabels = ["PBS", "DCZ"]
# XTickLabels = ["PBS", "PBS"]
# XTickLabels = ["PBS", "PBS", "DCZ"]
XLabel      = ""


#% Create figure and axis
Fig = plt.figure(figsize = [Fig_Width, Fig_Height], dpi = 300, label = mouse)

ax = plt.axes([Axis_left, Axis_bottom, Axis_width, Axis_height])   


# Single mice lines
SingleMice_Lines = ax.plot(Days2Use_int,
                           Data2Use.T,
                           c = SingleMice_Color,
                           marker = SingleMice_Marker,
                           markersize = SingleMice_MarkerSize,
                           linewidth = SingleMice_LineWidth)


# Mean error bars
MeanErrorBars_Graph = ax.errorbar(Days2Use_int,
                                np.mean(Data2Use, axis = 0),
                                c = MeanErrorBars_Color,
                                linewidth = MeanErrorBars_LineWidth,
                                yerr = np.std(Data2Use, axis = 0)/2,
                                ecolor = MeanErrorBars_Color,
                                elinewidth = MeanErrorBars_LineWidth,
                                capsize = 2.5,
                                capthick = MeanErrorBars_LineWidth,)




#Y axis
ax.set_ylabel(YLabel, fontsize = FontSize_AxisLabel, fontfamily = FontName)
ax.set_ylim(YLim)
ax.set_yticks(YTicks)
ax.set_yticklabels(YTickLabels, fontsize = FontSize_TickLabel, fontfamily = FontName)



ax.set_xlim(X_Lim)
ax.set_xticks(XTicks)
ax.set_xticklabels(XTickLabels, 
                   fontsize = FontSize_TickLabel, 
                   fontfamily = FontName,)

ax.set_xlabel(XLabel, fontsize = FontSize_AxisLabel, fontfamily = FontName)

# Remove spines
ax.spines[['top', "right"]].set_visible(False)

ax.spines[['left', "bottom"]].set_linewidth(Axis_LineWidth)


# increase tick width
ax.tick_params(width = Axis_LineWidth)

    
    
# legend
if Legend_ShowFlag:
    ax.legend(SingleMice_Lines, 
              ["A", "B", "C", "D"],
              # bbox_to_anchor = [0, .9, .9, .01],
              loc = [Legend_X0, Legend_Y0], 
              ncols = 2,
              fontsize = FontSize_TickLabel,
              prop = {"family" : FontName},
              frameon = False,
              fancybox = False,
              edgecolor = None,
              columnspacing = Legend_ColumnSpacing,
              handlelength = Legend_handlelength,
              handleheight = .2)

# Save
if Fig_SaveFlag:
    
    CreateDestinationFolder(Destination_Path = Fig_SaveDir)
    
    plt.savefig(os.path.join(Fig_SaveDir, Fig_SaveName), dpi='figure')
    
    
    
# %%Test significance
Contrast2Plot = 5
print("\nContrast %s%%" % (Contrast2Plot))


# ttest_ind  ttest_rel
Res = stats.ttest_rel(np.mean(Data2Use_PBS, axis = 1),
                      np.mean(Data2Use_DCZ, axis = 1),
                      alternative = "greater")

print("\npValue Mean Days PBS Greater than DCZ: %.4f" % (Res.pvalue))# greater less  two-sided
    
# 
Res = stats.ttest_rel(np.mean(Data2Use_DCZ, axis = 1),
                      np.mean(Data2Use_PBS, axis = 1),
                      alternative = "greater")


print("\npValue Mean days DCZ Greater than PBS: %.4f" % (Res.pvalue))# greater less  two-sided
    

# Mean and STD

print("\n\nPBS Percentage dot interceptions; mean +- std:\n%.2f +- %.2f" % (np.mean(np.mean(Data2Use_PBS, axis = 1)) * 100, 
                                                                            np.std(np.mean(Data2Use_PBS, axis = 1)) * 100))


print("\nDCZ Percentage dot interceptions; mean +- std:\n%.2f +- %.2f" % (np.mean(np.mean(Data2Use_DCZ, axis = 1)) * 100, 
                                                                            np.std(np.mean(Data2Use_DCZ, axis = 1)) * 100))
