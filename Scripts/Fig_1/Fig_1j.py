# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 13:53:51 2024

@author: Daniel
"""



import os
import pickle
import numpy as np
import matplotlib.pyplot as plt

# plt.close("all")
#%%
# DB_dir  = r"O:\OneDrive - MRC Laboratory of Molecular Biology\Projects\Ongoing\Inhibition\Imaging\VGAT\Data\DB"
DB_dir  = r"O:\OneDrive - MRC Laboratory of Molecular Biology\Work\Projects\Ongoing\Inhibition\Experiments\Completed\Imaging\VGAT\Data\DB"
# DB_name = "DB_VGAT_Dark_Smooth_2"
DB_name = "DB_VGAT_SingleRecs_Dark_Smooth"


with open(os.path.join(DB_dir, DB_name), 'rb') as handle:
    DB_Recs = pickle.load(handle)
    
    
Recs = list(DB_Recs)
nRecs = len(Recs)


#%% Get pearson correlation coefficients
Thres_ZScore = 5
Thres_Ratio = 0.5   # 0.5#.5



Pearson_Gyro_Tonic   = np.empty(0)
Pearson_Gyro_Phasic  = np.empty(0)
for rec in Recs:
    
    # Select responding neurons
    Responding_Neurons = np.flatnonzero(DB_Recs[rec]["ZScore_Max"] >= Thres_ZScore)
    
    # Get DF traces
    X = DB_Recs[rec]["DFoverF"][Responding_Neurons, :].T    
    
    nTimePoints, nNeurons = X.shape
    
    
    # Get middle of traces range
    Middle = (np.max(X, axis = 0) + np.min(X, axis = 0) ) / 2
    
    
    # Divide timepoints into active and inactive
    TimePoints_Neurons_Active   = X > Middle
    TimePoints_Neurons_Inactive = X < Middle
    
    
    # Compute ratio time spent active on inactive
    Ratio_Active_Inactive = np.sum(TimePoints_Neurons_Active, axis = 0) / np.sum(TimePoints_Neurons_Inactive, axis = 0)
    
    
    # Classify neurons into tonic or phasic
    TonicNeurons    = np.flatnonzero(Ratio_Active_Inactive > Thres_Ratio)
    PhasicNeurons   = np.flatnonzero(Ratio_Active_Inactive < Thres_Ratio)

    
    # Tonic
    for neuron in TonicNeurons:
        
        Pearson_Gyro_Tonic = np.append(Pearson_Gyro_Tonic, 
                                       np.corrcoef(X[:, neuron], 
                                                   np.abs(DB_Recs[rec]["Gyro_Smooth"]), 
                                                   rowvar = False
                                                   )[0, 1]
                                       )
        
    
    # Phasic
    for neuron in PhasicNeurons:
        
        Pearson_Gyro_Phasic = np.append(Pearson_Gyro_Phasic, 
                                        np.corrcoef(X[:, neuron], 
                                                    np.abs(DB_Recs[rec]["Gyro_Smooth"]), 
                                                    rowvar = False
                                                    )[0, 1]
                                        )
                                        

print("Tonic Neuron: Median: %.2f, Interquartile range: [%.2f, %.2f]" % (np.percentile(Pearson_Gyro_Tonic, 50), 
                                                                       np.percentile(Pearson_Gyro_Tonic, 25),
                                                                       np.percentile(Pearson_Gyro_Tonic, 75)))

print("Phasic Neuron: Median: %.2f, Interquartile range: [%.2f, %.2f]" % (np.percentile(Pearson_Gyro_Phasic, 50), 
                                                                          np.percentile(Pearson_Gyro_Phasic, 25),
                                                                          np.percentile(Pearson_Gyro_Phasic, 75)))


#%% Box plot Coeff Pearson Abs Gyro

save_flag   = False
save_dir    = r"O:\OneDrive - MRC Laboratory of Molecular Biology\General - VGAT Paper\Figures\Paper\Figure1\CorrelationTonic_Phasic"
save_name   = "BoxPlot_CorrCoeff_Tonic_Phasic_gyro.svg"


Color_Tonic         = (0, 0.68, 0.94)
Color_Phasic        = (0.93, 0, 0.55)

Colors = [Color_Tonic, Color_Phasic]


#Box plot params
Notch_Flag      = True
bootstrap       = 10000
whis            = (25, 75) # (5, 95)
Box_Positions   = [0, .5]
space           = .2


# X and Y axes tick and label
# XLabel      = "Comparaison between neurons"
XLabel      = "neurons"
Xticklabels = ["Tonic", "Phasic"]

YLim        = [np.percentile(Pearson_Gyro_Tonic, 24), 
               np.percentile(Pearson_Gyro_Phasic, 76)]
YTicks      = [-0.25, 0, .25]
YLabel      = "Pearson Coeff. with Gyro"


# Figure
Fig_Width   = 1
Fig_Height  = 1.2


# Axis
Axis_left       = .44
Axis_bottom     = 0.22
Axis_width      = .55
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
                        Pearson_Gyro_Tonic, 
                        Pearson_Gyro_Phasic,                      
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
xlim = [Box_Positions[0] - space, Box_Positions[1] + space]
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
