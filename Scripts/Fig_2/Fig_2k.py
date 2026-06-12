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
DF_Dir = r"O:\OneDrive - MRC Laboratory of Molecular Biology\Work\Projects\Ongoing\Inhibition\Experiments\Completed\Imaging\All\Code\MB_GR"

df_All = pd.read_pickle(os.path.join(DF_Dir, "df_All"))


#%% MotorTuned cells
Thres_PValue        = .05
Thres_PearsonCoeff  = 0.1 # .1
Thres_ZScore        = 5
Thres_QI            = 0.3


df_MotorTuned = df_All.loc[ (df_All["Gyro_Dark"]["ZScore_Max"] >= Thres_ZScore) *\
                            (df_All["Gyro_Dark"]["PearsonCoeff_absGyro"].abs() >= Thres_PearsonCoeff) *\
                            (df_All["Gyro_Dark"]["PearsonCoeff_absGyro_PValue"] <= Thres_PValue)
                            ]
    
df_MotorTuned_VisualResponsive = df_MotorTuned.loc[(df_MotorTuned["GR"]["QI"] >= Thres_QI) + \
                                                   (df_MotorTuned["MB"]["QI"] >= Thres_QI)
                                                   ]
    
    
#%% visual responsive motor tuned neurons only
# # SST
# Min_DFoverF_SST = pd.concat([df_MotorTuned_VisualResponsive.loc["SST", ("GR","Min_DFoverF")],
#                              df_MotorTuned_VisualResponsive.loc["SST", ("MB","Min_DFoverF")]], 
#                            axis = 1)    

# Min_DFoverF_SST = np.nanmin(Min_DFoverF_SST.values, axis = 1)


    
# # CCK
# Min_DFoverF_CCK = pd.concat([df_MotorTuned_VisualResponsive.loc["CCK", ("GR","Min_DFoverF")],
#                              df_MotorTuned_VisualResponsive.loc["CCK", ("MB","Min_DFoverF")]], 
#                            axis = 1)

# Min_DFoverF_CCK = np.nanmin(Min_DFoverF_CCK.values, axis = 1)


# # PV
# Min_DFoverF_PV = pd.concat([df_MotorTuned_VisualResponsive.loc["PV", ("GR","Min_DFoverF")],
#                             df_MotorTuned_VisualResponsive.loc["PV", ("MB","Min_DFoverF")]], 
#                            axis = 1)

# Min_DFoverF_PV = np.nanmin(Min_DFoverF_PV.values, axis = 1)


#%% All motor tuned
# SST
Min_DFoverF_SST = pd.concat([df_MotorTuned.loc["SST", ("GR","Min_DFoverF")],
                             df_MotorTuned.loc["SST", ("MB","Min_DFoverF")]], 
                            axis = 1)    

Min_DFoverF_SST = np.nanmin(Min_DFoverF_SST.values, axis = 1)


    
# CCK
Min_DFoverF_CCK = pd.concat([df_MotorTuned.loc["CCK", ("GR","Min_DFoverF")],
                             df_MotorTuned.loc["CCK", ("MB","Min_DFoverF")]], 
                            axis = 1)

Min_DFoverF_CCK = np.nanmin(Min_DFoverF_CCK.values, axis = 1)


# PV
Min_DFoverF_PV = pd.concat([df_MotorTuned.loc["PV", ("GR","Min_DFoverF")],
                            df_MotorTuned.loc["PV", ("MB","Min_DFoverF")]], 
                           axis = 1)

Min_DFoverF_PV = np.nanmin(Min_DFoverF_PV.values, axis = 1)

# Remove NAN
Min_DFoverF_PV = Min_DFoverF_PV[~np.isnan(Min_DFoverF_PV)]


#%% ‘two-sided’ two-sample Kolmogorov-Smirnov test
# two-sided: The null hypothesis is that the two distributions are identical, F(x)=G(x) for all x; 
# the alternative is that they are not identical. 
# The statistic is the maximum absolute difference between the empirical distribution functions of the samples.

Res_PV_SST = scipy.stats.ks_2samp(Min_DFoverF_PV, 
                                  Min_DFoverF_SST, 
                                  alternative = 'less')

print("\npValue two-sample KS test 1 tail PV less than SST:", Res_PV_SST.pvalue)


Res_PV_CCK = scipy.stats.ks_2samp(Min_DFoverF_PV, 
                                  Min_DFoverF_CCK, 
                                  alternative = 'two-sided')

print("\npValue two-sample KS test 2 tail PV different than CCK:", Res_PV_CCK.pvalue)


Res_SST_CCK = scipy.stats.ks_2samp(Min_DFoverF_CCK, 
                                   Min_DFoverF_SST, 
                                   alternative = 'less')

print("\npValue two-sample KS test 1 tail CCK less than SST:", Res_SST_CCK.pvalue)


#%% LogScale
save_fig    = False
save_dir    = r"O:\OneDrive - MRC Laboratory of Molecular Biology\General - VGAT Paper\Figures\Paper\Figure4\Ver5\CDFs"
save_name   = "CDF_MaxNegResp_MotorTuned.svg"



XMin = np.nanmin(np.hstack((Min_DFoverF_PV, Min_DFoverF_SST, Min_DFoverF_CCK)))
XMax = np.nanmax(np.hstack((Min_DFoverF_PV, Min_DFoverF_SST, Min_DFoverF_CCK)))


# Param figure
BinSize             = .5
rwidth              = .7
alpha               = .6
Traces_linewidth    = .8
line_linewidth      = .6


Color_SST   = "#7570b3"
Color_CCK   = "#d95f02"
Color_PV    = "#1b9e77"
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


n_SST, _, patches = ax.hist(Min_DFoverF_SST,
                        bins = Bin_Edges, 
                        density = True, 
                        cumulative = True,
                        histtype = 'step',
                        color = Color_SST)

n_CCK, _, patches = ax.hist(Min_DFoverF_CCK,
                        bins = Bin_Edges, 
                        density = True, 
                        cumulative = True,
                        histtype = 'step',
                        color = Color_CCK)

n_PV, _, patches = ax.hist(Min_DFoverF_PV,
                        bins = Bin_Edges, 
                        density = True, 
                        cumulative = True,
                        histtype = 'step',
                        color = Color_PV)



Fig = plt.figure(figsize = [Fig_Width, Fig_Height], dpi = 300)

ax = plt.axes([Axis_left, Axis_bottom, Axis_width, Axis_height])   


ax.plot(Bin_Centres, n_SST, linewidth = Traces_linewidth, c = Color_SST)
ax.plot(Bin_Centres, n_CCK, linewidth = Traces_linewidth, c = Color_CCK)
ax.plot(Bin_Centres, n_PV, linewidth = Traces_linewidth, c = Color_PV)


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
    
    