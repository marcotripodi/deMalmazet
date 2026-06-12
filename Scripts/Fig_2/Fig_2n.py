# -*- coding: utf-8 -*-
"""
Created on Sun Mar  8 21:34:55 2026

@author: Daniel
"""


import os
import numpy as np
import matplotlib.pyplot as plt
from PathGeneral_Func import CreateDestinationFolder
import pandas as pd
from tqdm import tqdm


#%%
DF_Dir = r"O:\OneDrive - MRC Laboratory of Molecular Biology\Work\Projects\Ongoing\Inhibition\Experiments\Completed\Imaging\All\Code\MB\MovingBarsRF"
            
# df_All = pd.read_pickle(os.path.join(DF_Dir, "df_All"))
df_All_Intersect = pd.read_pickle(os.path.join(DF_Dir, "df_All_Intersect"))


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
df_MotorTuned = df_All_Intersect.loc[ (df_All_Intersect["Gyro_Dark"]["ZScore_Max"] >= Thres_ZScore) *\
                                     (df_All_Intersect["Gyro_Dark"]["PearsonCoeff_absGyro"].abs() >= Thres_PearsonCoeff) *\
                                         (df_All_Intersect["Gyro_Dark"]["PearsonCoeff_absGyro_PValue"] <= Thres_PValue)
                                         ]
    
df_MotorTuned_SST = df_MotorTuned.loc["SST"]

df_MotorTuned_CCK = df_MotorTuned.loc["CCK"]

df_MotorTuned_PV = df_MotorTuned.loc["PV"]


#%%
RFSizes_MotorTuned_SST = df_MotorTuned_SST["MB", "RFsize_deg"][df_MotorTuned_SST["MB", "QI"] > Thres_QI]
# df_MotorTuned_SST.loc[slice(None), ("MB", "QI")]

RFSizes_MotorTuned_CCK = df_MotorTuned_CCK["MB", "RFsize_deg"][df_MotorTuned_CCK["MB", "QI"] > Thres_QI]

RFSizes_MotorTuned_PV = df_MotorTuned_PV["MB", "RFsize_deg"][df_MotorTuned_PV["MB", "QI"] > Thres_QI]


#%% Stats

print("RF size SST: Median: %.2f, Interquartile range: [%.2f, %.2f]" % (np.percentile(RFSizes_MotorTuned_SST, 50), 
                                                                       np.percentile(RFSizes_MotorTuned_SST, 25),
                                                                       np.percentile(RFSizes_MotorTuned_SST, 75)))

print("# n neurons SST: %i\n" % (len(RFSizes_MotorTuned_SST)))


print("RF size CCK: Median: %.2f, Interquartile range: [%.2f, %.2f]" % (np.percentile(RFSizes_MotorTuned_CCK, 50), 
                                                                       np.percentile(RFSizes_MotorTuned_CCK, 25),
                                                                       np.percentile(RFSizes_MotorTuned_CCK, 75)))

print("# n neurons CCK: %i\n" % (len(RFSizes_MotorTuned_CCK)))


print("RF size PV: Median: %.2f, Interquartile range: [%.2f, %.2f]" % (np.percentile(RFSizes_MotorTuned_PV, 50), 
                                                                       np.percentile(RFSizes_MotorTuned_PV, 25),
                                                                       np.percentile(RFSizes_MotorTuned_PV, 75)))

print("# n neurons PV: %i\n" % (len(RFSizes_MotorTuned_PV)))


#%% Box plot Coeff Pearson Abs Gyro
save_flag   = False
save_dir    = r"O:\OneDrive - MRC Laboratory of Molecular Biology\General - VGAT Paper\Figures\Paper\Figure4\Ver5\RF"
save_name   = "BoxPlot_RFSizes.svg"


Color_SST   = "#7570b3"
Color_CCK   = "#d95f02"
Color_PV    = "#1b9e77"

Colors = [Color_SST, Color_CCK, Color_PV]


#Box plot params
Notch_Flag      = False
bootstrap       = 10000
whis            = (25, 75) # (5, 95)
Box_Positions   = [0, .5, 1]
space           = .2


# X and Y axes tick and label
XLabel      = ""
Xticklabels = ["SST\nVGAT", "CCK\nVGAT", "PV\nVGAT"]

YLim        = [0, 31]
YTicks      = [0, 15, 30]
YLabel      = "Receptive Field size (\u00b0)"


# Figure
Fig_Width   = 1.3
Fig_Height  = 1.1


# Axis
Axis_left       = .35
Axis_bottom     = 0.18
Axis_width      = .63
Axis_height     = .65
Axis_LineWidth  = .5


FontSize_TickLabel  = 6
FontSize_AxisLabel  = 7
FontName            = "arial"


# Create figure
Fig = plt.figure(figsize = [Fig_Width, Fig_Height], dpi = 300)

# Create axis
ax = plt.axes([Axis_left, Axis_bottom, Axis_width, Axis_height])

# Plot box plot
bplot = ax.boxplot( x = [
                            RFSizes_MotorTuned_SST, 
                            RFSizes_MotorTuned_CCK, 
                            RFSizes_MotorTuned_PV
                            ],
                        sym = "",
                        patch_artist = True,
                        notch = Notch_Flag, 
                        bootstrap = bootstrap, 
                        whis = whis,
                        positions = Box_Positions
                        )


# Draw Significance line
# ax.plot(ax.get_xlim(), [PVal_thres, PVal_thres], linewidth=PVal_thres_linewidth, c=PVal_thres_Color)


# Adjust xticks and labels
xlim = [Box_Positions[0] - space, Box_Positions[-1] + space]
ax.set_xlim(xlim)
ax.set_xticklabels(Xticklabels, fontsize = FontSize_TickLabel, fontfamily = FontName)
ax.set_xlabel(XLabel, fontsize = FontSize_AxisLabel, fontfamily = FontName)
ax.xaxis.set_tick_params(length=0, width=0, which='major')

ax.set_ylim(YLim)
ax.set_yticks(YTicks)
ax.set_yticklabels(YTicks, fontsize=FontSize_TickLabel, fontfamily=FontName)
ax.set_ylabel(YLabel, fontsize=FontSize_AxisLabel, fontfamily=FontName)


# Remove spines
ax.spines[['top', "right", 'bottom']].set_visible(False)

ax.spines[['left']].set_linewidth(Axis_LineWidth)


# increase tick width
ax.tick_params(width = Axis_LineWidth)


# fill with colors
for patch, color in zip(bplot['boxes'], Colors):
    patch.set_facecolor(color)
    # patch.set_edgecolor(color)
    patch.set_linewidth(.5)
    

# for patch, color in zip(bplot['whiskers'], Colors):
for patch in bplot['whiskers']:
    patch.set_color("k")
    patch.set_linewidth(.5)
    

# for patch, color in zip(bplot['caps'], Colors):
for patch in bplot['caps']:
    patch.set_color("k")
    patch.set_linewidth(.5)
    
    
for patch, color in zip(bplot['medians'], Colors):
    patch.set_color("k")
    patch.set_linewidth(.5)
    
    
# save
if save_flag:
    plt.savefig(os.path.join(save_dir, save_name), dpi='figure')
