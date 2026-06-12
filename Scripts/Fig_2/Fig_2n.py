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
DF_Dir = r""
            
df_All_Intersect = pd.read_pickle(os.path.join(DF_Dir, "df_All_Intersect"))


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
df_MotorTuned = df_All_Intersect.loc[ (df_All_Intersect["Gyro_Dark"]["ZScore_Max"] >= Thres_ZScore) *\
                                     (df_All_Intersect["Gyro_Dark"]["PearsonCoeff_absGyro"].abs() >= Thres_PearsonCoeff) *\
                                         (df_All_Intersect["Gyro_Dark"]["PearsonCoeff_absGyro_3alue"] <= Thres_3alue)
                                         ]
    
df_MotorTuned_1 = df_MotorTuned.loc["1"]

df_MotorTuned_2 = df_MotorTuned.loc["2"]

df_MotorTuned_3 = df_MotorTuned.loc["3"]


#%%
RFSizes_MotorTuned_1 = df_MotorTuned_1["MB", "RFsize_deg"][df_MotorTuned_1["MB", "QI"] > Thres_QI]
# df_MotorTuned_1.loc[slice(None), ("MB", "QI")]

RFSizes_MotorTuned_2 = df_MotorTuned_2["MB", "RFsize_deg"][df_MotorTuned_2["MB", "QI"] > Thres_QI]

RFSizes_MotorTuned_3 = df_MotorTuned_3["MB", "RFsize_deg"][df_MotorTuned_3["MB", "QI"] > Thres_QI]


#%% Stats

print("RF size 1: Median: %.2f, Interquartile range: [%.2f, %.2f]" % (np.percentile(RFSizes_MotorTuned_1, 50), 
                                                                       np.percentile(RFSizes_MotorTuned_1, 25),
                                                                       np.percentile(RFSizes_MotorTuned_1, 75)))

print("# n neurons 1: %i\n" % (len(RFSizes_MotorTuned_1)))


print("RF size 2: Median: %.2f, Interquartile range: [%.2f, %.2f]" % (np.percentile(RFSizes_MotorTuned_2, 50), 
                                                                       np.percentile(RFSizes_MotorTuned_2, 25),
                                                                       np.percentile(RFSizes_MotorTuned_2, 75)))

print("# n neurons 2: %i\n" % (len(RFSizes_MotorTuned_2)))


print("RF size 3: Median: %.2f, Interquartile range: [%.2f, %.2f]" % (np.percentile(RFSizes_MotorTuned_3, 50), 
                                                                       np.percentile(RFSizes_MotorTuned_3, 25),
                                                                       np.percentile(RFSizes_MotorTuned_3, 75)))

print("# n neurons 3: %i\n" % (len(RFSizes_MotorTuned_3)))


#%% Box plot Coeff Pearson Abs Gyro
save_flag   = False
save_dir    = r""
save_name   = "BoxPlot_RFSizes.svg"


Color_1   = "#7570b3"
Color_2   = "#d95f02"
Color_3    = "#1b9e77"

Colors = [Color_1, Color_2, Color_3]


#Box plot params
Notch_Flag      = False
bootstrap       = 10000
whis            = (25, 75) # (5, 95)
Box_Positions   = [0, .5, 1]
space           = .2


# X and Y axes tick and label
XLabel      = ""
Xticklabels = ["1\nVGAT", "2\nVGAT", "3\nVGAT"]

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
                            RFSizes_MotorTuned_1, 
                            RFSizes_MotorTuned_2, 
                            RFSizes_MotorTuned_3
                            ],
                        sym = "",
                        patch_artist = True,
                        notch = Notch_Flag, 
                        bootstrap = bootstrap, 
                        whis = whis,
                        positions = Box_Positions
                        )


# Draw Significance line
# ax.plot(ax.get_xlim(), [3al_thres, 3al_thres], linewidth=3al_thres_linewidth, c=3al_thres_Color)


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
