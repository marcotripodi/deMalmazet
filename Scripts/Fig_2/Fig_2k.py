# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 13:53:51 2024

@author: Daniel
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from PathGeneral_Func import CreateDestinationFolder
import scipy


#%%
DF_Dir = r""

df_All = pd.read_pickle(os.path.join(DF_Dir, "df_All"))


#%% MotorTuned cells
Thres_Pvalue        = .05
Thres_PearsonCoeff  = 0.1 # .1
Thres_ZScore        = 5
Thres_QI            = 0.3


df_MotorTuned = df_All.loc[ (df_All["Gyro_Dark"]["ZScore_Max"] >= Thres_ZScore) *\
                            (df_All["Gyro_Dark"]["PearsonCoeff_absGyro"].abs() >= Thres_PearsonCoeff) *\
                            (df_All["Gyro_Dark"]["PearsonCoeff_absGyro_3alue"] <= Thres_Pvalue)
                            ]
    
df_MotorTuned_VisualResponsive = df_MotorTuned.loc[(df_MotorTuned["GR"]["QI"] >= Thres_QI) + \
                                                   (df_MotorTuned["MB"]["QI"] >= Thres_QI)
                                                   ]


#%% All motor tuned
# 1
Min_DFoverF_1 = pd.concat([df_MotorTuned.loc["Marker1", ("GR","Min_DFoverF")],
                             df_MotorTuned.loc["Marker1", ("MB","Min_DFoverF")]], 
                            axis = 1)    

Min_DFoverF_1 = np.nanmin(Min_DFoverF_1.values, axis = 1)


    
# 2
Min_DFoverF_2 = pd.concat([df_MotorTuned.loc["2", ("GR","Min_DFoverF")],
                             df_MotorTuned.loc["2", ("MB","Min_DFoverF")]], 
                            axis = 1)

Min_DFoverF_2 = np.nanmin(Min_DFoverF_2.values, axis = 1)


# 3
Min_DFoverF_3 = pd.concat([df_MotorTuned.loc["3", ("GR","Min_DFoverF")],
                            df_MotorTuned.loc["3", ("MB","Min_DFoverF")]], 
                           axis = 1)

Min_DFoverF_3 = np.nanmin(Min_DFoverF_3.values, axis = 1)

# Remove NAN
Min_DFoverF_3 = Min_DFoverF_3[~np.isnan(Min_DFoverF_3)]


#%% ‘two-sided’ two-sample Kolmogorov-Smirnov test
# two-sided: The null hypothesis is that the two distributions are identical, F(x)=G(x) for all x; 
# the alternative is that they are not identical. 
# The statistic is the maximum absolute difference between the empirical distribution functions of the samples.

Res_3_1 = scipy.stats.ks_2samp(Min_DFoverF_3, 
                                  Min_DFoverF_1, 
                                  alternative = 'less')

print("\nPvalue two-sample KS test 1 tail 3 less than 1:", Res_3_1.Pvalue)


Res_3_2 = scipy.stats.ks_2samp(Min_DFoverF_3, 
                                  Min_DFoverF_2, 
                                  alternative = 'two-sided')

print("\nPvalue two-sample KS test 2 tail 3 different than 2:", Res_3_2.Pvalue)


Res_1_2 = scipy.stats.ks_2samp(Min_DFoverF_2, 
                                   Min_DFoverF_1, 
                                   alternative = 'less')

print("\nPvalue two-sample KS test 1 tail 2 less than 1:", Res_1_2.Pvalue)


#%% LogScale
save_fig    = False
save_dir    = r""
save_name   = ".svg"



XMin = np.nanmin(np.hstack((Min_DFoverF_3, Min_DFoverF_1, Min_DFoverF_2)))
XMax = np.nanmax(np.hstack((Min_DFoverF_3, Min_DFoverF_1, Min_DFoverF_2)))


# Param figure
BinSize             = .5
rwidth              = .7
alpha               = .6
Traces_linewidth    = .8
line_linewidth      = .6


Color_1   = "#7570b3"
Color_2   = "#d95f02"
Color_3    = "#1b9e77"
Color_line  = ".65"


FontSize_TickLabel  = 6
FontSize_AxisLabel  = 7
FontName            = "arial"


XTicks = [-1, -.5, 0]
XLabel = "Min DF/F MB GR"

YTicks  = [0, .5, 1]
YLabel  = "Fraction of \nmotor tuned neurons"
# YLabel  = ""
YLim    = [0, 1.01]

# figure
Fig_Width = 1.5
Fig_Height = 1.2


# Axis
Axis_left       = .25
Axis_bottom     = 0.3
Axis_width      = .7
Axis_height     = .6
Axis_LineWidth  = .5


Bins_num = 30
Bin_Edges = np.linspace(XMin, XMax, num = Bins_num)

Bin_Centres = (Bin_Edges[:-1] + Bin_Edges[1:])/2
#

#% Plot the PDF.
Fig = plt.figure(figsize = [Fig_Width, Fig_Height], dpi = 300)

ax = plt.axes([Axis_left, Axis_bottom, Axis_width, Axis_height])   


n_1, _, patches = ax.hist(Min_DFoverF_1,
                        bins = Bin_Edges, 
                        density = True, 
                        cumulative = True,
                        histtype = 'step',
                        color = Color_1)

n_2, _, patches = ax.hist(Min_DFoverF_2,
                        bins = Bin_Edges, 
                        density = True, 
                        cumulative = True,
                        histtype = 'step',
                        color = Color_2)

n_3, _, patches = ax.hist(Min_DFoverF_3,
                        bins = Bin_Edges, 
                        density = True, 
                        cumulative = True,
                        histtype = 'step',
                        color = Color_3)



Fig = plt.figure(figsize = [Fig_Width, Fig_Height], dpi = 300)

ax = plt.axes([Axis_left, Axis_bottom, Axis_width, Axis_height])   


ax.plot(Bin_Centres, n_1, linewidth = Traces_linewidth, c = Color_1)
ax.plot(Bin_Centres, n_2, linewidth = Traces_linewidth, c = Color_2)
ax.plot(Bin_Centres, n_3, linewidth = Traces_linewidth, c = Color_3)


# Plot .5 line
# ax.plot([Thres_QI, Thres_QI], [0, 1], 
#         linewidth = line_linewidth, 
#         linestyle = "--", 
#         c = Color_line)


# x axis  
# ax.set_xscale('log', subs = [0])

ax.set_xlabel(XLabel, fontsize = FontSize_AxisLabel, fontfamily = FontName)
ax.set_xlim([Bin_Edges[0], Bin_Edges[-1]])
ax.set_xticks(XTicks)
ax.set_xticklabels(XTicks, fontsize = FontSize_TickLabel, fontfamily = FontName)


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
    
    