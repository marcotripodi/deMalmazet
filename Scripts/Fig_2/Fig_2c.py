# -*- coding: utf-8 -*-
"""
Created on Tue Dec 31 16:54:24 2024

@author: Daniel
"""

import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt


#%% Load data
Data_dir = r"O:\OneDrive - MRC Laboratory of Molecular Biology\Work\Projects\Ongoing\Inhibition\Experiments\Completed\RNAscope_Antonio\Code\Data"
Data_Name = "Overlap_DF"


Data_Path = os.path.join(Data_dir, Data_Name)

df = pd.read_pickle(Data_Path)



#%% Function to get percentage mean and std for each population of overlap
def GetPercentageLayer(MainVGATPop, SubVGATPop, Label):
    
    """
    SubVGATPop = "# PV x SST x VGAT"
    
    MainVGATPop = '# PV x VGAT'
    
    Label = "% SST in PVxVGAT"
    """
    
    
    Intersect = df.loc[df[SubVGATPop].notna(), 
                     ['layer', MainVGATPop, SubVGATPop]]


    Percentage = pd.DataFrame({"layer" : Intersect["layer"], 
                               Label: Intersect[SubVGATPop] / Intersect[MainVGATPop] * 100}
                 )
    
    
    Mean_STD = np.empty((3, 2))
    i = 0
    for layer in ["SCs", "SCm", "SCd"]:
        
        Mean_STD[i, 0] =  Percentage.loc[Percentage["layer"] == layer, 
                                         Label].mean()
        
                
        Mean_STD[i, 1] =  Percentage.loc[Percentage["layer"] == layer,
                                         Label].std()
    
        i += 1
        
        
    return pd.DataFrame(data = Mean_STD, 
                        index = ['SCs', 'SCm', 'SCd'], 
                        columns = [[Label, Label], ["mean", "std"]])


#%%
# CCK in SST VGAT
SubVGATPop  = '# CCK x SST x VGAT'
MainVGATPop = '# SST x VGAT'
Label       = "% CCK in SST x VGAT"
df_CCKinSST  = GetPercentageLayer(MainVGATPop, SubVGATPop, Label)
print(df_CCKinSST)


# PV in SST VGAT
SubVGATPop  = "# PV x SST x VGAT"
MainVGATPop = '# SST x VGAT'
Label       = "% PV in SST x VGAT"
df_PVinSST  = GetPercentageLayer(MainVGATPop, SubVGATPop, Label)
print(df_PVinSST)


#%% aligned vertical
save_fig    = False
save_dir    = r"O:\OneDrive - MRC Laboratory of Molecular Biology\General - VGAT Paper\Figures\Paper\Figure2\Ver7"
save_name   = "SST_VGAT_vert.svg"


# Param figure
BarWidth = .8
# Error Bars
capsize     = 3.5
capthick    = .4
elinewidth  = .4


Color_SST   = "#7570b3"
Color_CCK   = "#d95f02"
Color_PV    = "#1b9e77"


FontSize_TickLabel = 6
FontSize_AxisLabel = 7
FontName = "arial"


XTicks = []
XLabel = ""
YTicks = [0, 25, 50]
# YLabel = "% of $SST^{ON}$x $VGAT^{ON}$" # $\mathregular{H^{GO}}$
YLabel = "% of $\mathregular{SST^{ON}}$x $\mathregular{VGAT^{ON}}$" # $\mathregular{H^{GO}}$
xlim_slack = .1

# YLim = [0, np.max(df_CCKinPV["% CCK in PV x VGAT"])]
# YLim = [0, 23]
YLim = [0, 59]

# figure
Fig_Width = .75
Fig_Height = 1


# Axis
Axis_left       = .47
Axis_bottom     = 0.05
Axis_right      = .95
Axis_top        = .9
Axis_LineWidth  = .5


#% Plot the PDF.
fig, axs = plt.subplots(3, 1,
                        figsize = [Fig_Width, Fig_Height], 
                        dpi = 300,
                        sharex = True,
                        sharey = True,
                        gridspec_kw = {"left": Axis_left, 
                                       "right": Axis_right, 
                                       "top": Axis_top, 
                                       "bottom": Axis_bottom}
                        )

x = [0, 1]

for ax, SClayer in zip(axs, df_CCKinSST.index):
    
    ax.bar(0,
           df_CCKinSST.loc[SClayer, ("% CCK in SST x VGAT", "mean")],
           width = BarWidth,
           color = Color_CCK,
           yerr = df_CCKinSST.loc[SClayer, ("% CCK in SST x VGAT", "std")],
           ecolor = Color_CCK,
           capsize = capsize,
           error_kw = {"elinewidth": elinewidth,
                       "capthick": capthick}
           )
    
    ax.bar(1,
           df_PVinSST.loc[SClayer, ("% PV in SST x VGAT", "mean")],
           width = BarWidth,
           color = Color_PV,
           yerr = df_PVinSST.loc[SClayer, ("% PV in SST x VGAT", "std")],
           ecolor = Color_PV,
           capsize = capsize,
           error_kw = {"elinewidth": elinewidth,
                       "capthick": capthick}
           )


    # x axis  
    ax.set_xlabel(XLabel, fontsize = FontSize_AxisLabel, fontfamily = FontName)
    ax.set_xlim([x[0] - BarWidth/2 - xlim_slack, 
                 x[-1] + BarWidth/2 + xlim_slack])
    
    ax.set_xticks(XTicks)
    ax.set_xticklabels(XTicks, fontsize = FontSize_TickLabel, fontfamily = FontName)
    
    
    #Y axis
    ax.set_yticks(YTicks)
    # ax.set_yticklabels(["", ""], fontsize = FontSize_TickLabel, fontfamily = FontName)
    ax.tick_params(axis='y', which='major', labelsize=FontSize_TickLabel)
    # ax.set_ylim([0, np.max(df_CCKinPV["% CCK in PV x VGAT"])])
    
    # Remove spines
    ax.spines[['top', 'bottom', "right"]].set_visible(False)
    
    ax.spines[['left']].set_linewidth(Axis_LineWidth)
    
    
    # increase tick width
    ax.tick_params(width = Axis_LineWidth)


# Ylabel    
axs[1].set_ylabel(YLabel, fontsize = FontSize_AxisLabel, fontfamily = FontName)
axs[1].set_ylim(YLim)
axs[1].set_yticklabels(YTicks, fontsize = FontSize_TickLabel, fontfamily = FontName)

# Save
if save_fig:
    
    # CreateDestinationFolder(Destination_Path = save_dir)
    save_path = os.path.join(save_dir, save_name)
    plt.savefig(save_path, dpi='figure')    
