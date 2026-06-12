# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 13:53:51 2024

@author: Daniel
"""

import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
from PathGeneral_Func import CreateDestinationFolder
import scipy


#%% Load DB
DB_dir  = r"O:\OneDrive - MRC Laboratory of Molecular Biology\Work\Projects\Ongoing\Inhibition\Experiments\Completed\Imaging\All\Data\DB_Baseline_Dark_Gyro"

# SST
# DB_dir  = r"O:\OneDrive - MRC Laboratory of Molecular Biology\Projects\Ongoing\Inhibition\Experiments\Completed\Imaging\SST\Data\DB"
DB_name = "DB_VGAT_SST_SingleRecs_Dark_Smooth"


with open(os.path.join(DB_dir, DB_name), 'rb') as handle:
    DB_Recs_SST = pickle.load(handle)
    
    
# CCK
# DB_dir  = r"O:\OneDrive - MRC Laboratory of Molecular Biology\Projects\Ongoing\Inhibition\Experiments\Completed\Imaging\CCK\Data\DB"
DB_name = "DB_VGAT_CCK_SingleRecs_Dark_Smooth"


with open(os.path.join(DB_dir, DB_name), 'rb') as handle:
    DB_Recs_CCK = pickle.load(handle)
    


# PV
# DB_dir  = r"O:\OneDrive - MRC Laboratory of Molecular Biology\Projects\Ongoing\Inhibition\Experiments\Completed\Imaging\PV\Data\DB"
DB_name = "DB_VGAT_PV_SingleRecs_Dark_Smooth"


with open(os.path.join(DB_dir, DB_name), 'rb') as handle:
    DB_Recs_PV = pickle.load(handle)
    
    

#%%
Thres_ZScore = 5


# SST
Ratio_Active_Inactiv_SST = np.empty(0)
for rec in DB_Recs_SST:
    
    # Select responding neurons
    Responding_Neurons = np.flatnonzero(DB_Recs_SST[rec]["ZScore_Max"] >= Thres_ZScore)
    
    
    # Get DF traces
    X = DB_Recs_SST[rec]["DFoverF"][Responding_Neurons, :].T    
    
    
    # Get middle of traces range♠
    Middle = (np.max(X, axis = 0) + np.min(X, axis = 0) ) / 2
    
    
    # Divide timepoints into active and inactive
    TimePoints_Neurons_Active   = X > Middle
    TimePoints_Neurons_Inactive = X < Middle
    
    
    # Compute ratio time spent active on inactive
    Ratio_Active_Inactiv_SST = np.append(Ratio_Active_Inactiv_SST, 
                                      np.sum(TimePoints_Neurons_Active, axis = 0) / np.sum(TimePoints_Neurons_Inactive, axis = 0))
    

# CCK
Ratio_Active_Inactiv_CCK = np.empty(0)
for rec in DB_Recs_CCK:
    
    # Select responding neurons
    Responding_Neurons = np.flatnonzero(DB_Recs_CCK[rec]["ZScore_Max"] >= Thres_ZScore)
    
    
    # Get DF traces
    X = DB_Recs_CCK[rec]["DFoverF"][Responding_Neurons, :].T    
    
    
    # Get middle of traces range
    Middle = (np.max(X, axis = 0) + np.min(X, axis = 0) ) / 2
    
    
    # Divide timepoints into active and inactive
    TimePoints_Neurons_Active   = X > Middle
    TimePoints_Neurons_Inactive = X < Middle
    
    
    # Compute ratio time spent active on inactive
    Ratio_Active_Inactiv_CCK = np.append(Ratio_Active_Inactiv_CCK, 
                                      np.sum(TimePoints_Neurons_Active, axis = 0) / np.sum(TimePoints_Neurons_Inactive, axis = 0))
    

# PV
Ratio_Active_Inactiv_PV = np.empty(0)
for rec in DB_Recs_PV:
    
    # Select responding neurons
    Responding_Neurons = np.flatnonzero(DB_Recs_PV[rec]["ZScore_Max"] >= Thres_ZScore)
    
    
    # Get DF traces
    X = DB_Recs_PV[rec]["DFoverF"][Responding_Neurons, :].T    
    
    
    # Get middle of traces range
    Middle = (np.max(X, axis = 0) + np.min(X, axis = 0) ) / 2
    
    
    # Divide timepoints into active and inactive
    TimePoints_Neurons_Active   = X > Middle
    TimePoints_Neurons_Inactive = X < Middle
    
    
    # Compute ratio time spent active on inactive
    Ratio_Active_Inactiv_PV = np.append(Ratio_Active_Inactiv_PV, 
                                      np.sum(TimePoints_Neurons_Active, axis = 0) / np.sum(TimePoints_Neurons_Inactive, axis = 0))
    
    
#%% ‘two-sided’ two-sample Kolmogorov-Smirnov test
# two-sided: The null hypothesis is that the two distributions are identical, F(x)=G(x) for all x; 
# the alternative is that they are not identical. 
# The statistic is the maximum absolute difference between the empirical distribution functions of the samples.

Res_PV_SST = scipy.stats.ks_2samp(Ratio_Active_Inactiv_PV, 
                                  Ratio_Active_Inactiv_SST, 
                                  alternative = 'two-sided')

Res_PV_SST.pvalue


Res_PV_CCK = scipy.stats.ks_2samp(Ratio_Active_Inactiv_PV, 
                                  Ratio_Active_Inactiv_CCK, 
                                  alternative = 'two-sided')

Res_PV_CCK.pvalue


Res_SST_CCK = scipy.stats.ks_2samp(Ratio_Active_Inactiv_SST, 
                                   Ratio_Active_Inactiv_CCK, 
                                   alternative = 'two-sided')

Res_SST_CCK.pvalue


print("pValue two-sample KS test SST - CCK:", Res_SST_CCK.pvalue)
print("pValue two-sample KS test PV - CCK:", Res_PV_CCK.pvalue)
print("pValue two-sample KS test SST - PV:", Res_PV_SST.pvalue)


#%% LogScale
save_fig    = False
save_dir    = r"O:\OneDrive - MRC Laboratory of Molecular Biology\General - VGAT Paper\Figures\Paper\Figure3\Ver5\Histogram_Ratios"
save_name   = "CDF_ratios_Logs_Tonic_Phasic_All.svg"


Thres_Ratio = .5


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


XTicks = [.005, .05, .5, 5]
# XTicks = [.001, .01, .1, .5, 10]
XLabel = "Ratio time active / inactive"

YTicks  = [0, .5, 1]
YLabel  = "Fraction of neurons"
YLim    = [0, 1.01]

# figure
Fig_Width = 1.5
Fig_Height = 1.2


# Axis
Axis_left       = .27
Axis_bottom     = 0.3
Axis_width      = .7
Axis_height     = .6
Axis_LineWidth  = .5


Bins_num = 30
# Bins = np.geomspace(Ratio_Active_Inactive.min(), Ratio_Active_Inactive.max(), num = Bins_num)
# Bins = np.geomspace(0.0018, 8.71, num = Bins_num)
Bin_Edges = np.geomspace(0.001, 10, num = Bins_num)

Bin_Centres = (Bin_Edges[:-1] + Bin_Edges[1:])/2
#

#% Plot the PDF. put 
Fig = plt.figure(figsize = [Fig_Width, Fig_Height], dpi = 300)#, visible = False)

ax = plt.axes([Axis_left, Axis_bottom, Axis_width, Axis_height])   


n_SST, _, patches = ax.hist(Ratio_Active_Inactiv_SST,
                        bins = Bin_Edges, 
                        density = True, 
                        cumulative = True,
                        histtype = 'step',
                        color = Color_SST)

n_CCK, _, patches = ax.hist(Ratio_Active_Inactiv_CCK,
                        bins = Bin_Edges, 
                        density = True, 
                        cumulative = True,
                        histtype = 'step',
                        color = Color_CCK)

n_PV, _, patches = ax.hist(Ratio_Active_Inactiv_PV,
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
ax.plot([Thres_Ratio, Thres_Ratio], [0, 1], 
        linewidth = line_linewidth, 
        linestyle = "--", 
        c = Color_line)


# x axis  
ax.set_xscale('log', subs = [0])

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
    
    