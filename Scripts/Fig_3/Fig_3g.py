# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 11:45:15 2024

@author: Daniel
"""

import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from PathGeneral_Func import CreateDestinationFolder


#%%

Counts_PV   = [28.9, 55.2, 44.6]

Counts_SST  = [30, 29.4, 29.2]

Counts_CCK  = [25.8, 16.1, 16.7]


df1 = pd.DataFrame( data = { "Counts" : Counts_SST,
                             "Gene": "SST"})


df2 = pd.DataFrame( data = { "Counts" : Counts_PV,
                             "Gene": "PV"})


df3 = pd.DataFrame( data = { "Counts" : Counts_CCK,
                             "Gene": "CCK"})



df = pd.concat([df1, df2, df3], ignore_index = True)


print("\nPercentage SSTonxVGATon among VGATon presynaptic to Pitx2 neurons. Mean +- STD\n%.2f +- %.2f" % (np.mean(Counts_SST), np.std(Counts_SST)))


#%% Strip_Plot
Save_Dir = r"O:\OneDrive - MRC Laboratory of Molecular Biology\General - VGAT Paper\Figures\Paper\Figure5\Ver1\Elena"
# Save_Dir = r"O:\OneDrive - MRC Laboratory of Molecular Biology\Projects\Ongoing\Inhibitory neurons SC" + \
#             r"\\RNAScope_Elena\Figures"
Save_Fig = False

Save_Name = "StriPlot_percentage subpopulation in VGAT pitx2 Pre.svg"


MarkerSize = 2.5


# palette = sns.color_palette("colorblind", n_colors = 2)
palette = sns.color_palette("colorblind", n_colors = 1)

Color_SST   = [141/255, 160/255, 203/255]
Color_PV    = [252/255, 141/255, 98/255]
Color_CCK   = [102/255, 194/255, 165/255]


# Line_PValue
Line_PValue_Color       = [.6, .6, .6]
Line_PValue_LineWidth   = .6
Line_PValue_LineStyle   = "--"


# X axis
XTicks      = [0]
xticklabels = ["$\mathregular{SST^{ON}}$"]
XLabel      = ""
XLim        = [-.1, .1]


# Y axis
YLabel          = "% of Inhibitory input\nto $\mathregular{Pitx2^{ON}}$ neurons"
yticks          = [0, 15, 30]
yticks_label    = yticks #["0", "10"]
YLim            = [0, 32]


# figure
Fig_Width   = 1
Fig_Height  = 1.2

    
# Axis
Axis_LineWidth  = .3
Axis_left       = .5
Axis_bottom     = 0.2
Axis_width      = .4
Axis_height     = .7


# Traces    
Linewidth_Traces    = .8
 

FontSize_TickLabel  = 6
FontSize_AxisLabel  = 7
FontName            = "arial"



Fig = plt.figure(figsize = [Fig_Width, Fig_Height], dpi = 300)

ax = plt.axes([Axis_left, Axis_bottom, Axis_width, Axis_height])   
    

ax = sns.swarmplot( ax = ax, 
                    data = df1,
                    x = "Gene",
                    y = "Counts",
                    hue = None,
                    color = "k",
                    size = MarkerSize,
                    linewidth = Axis_LineWidth,
                    legend = False, 
                    edgecolor = "none")


# Y axis
ax.set_ylabel(YLabel, fontsize = FontSize_AxisLabel, fontfamily = FontName)
ax.set_ylim(YLim)
ax.set_yticks(yticks)
ax.set_yticklabels(yticks_label, fontsize=FontSize_TickLabel, fontfamily=FontName)


# X axis
ax.set_xlim(XLim)
ax.set_xticks(XTicks)
ax.set_xticklabels(xticklabels, fontsize=FontSize_TickLabel, fontfamily=FontName)
ax.set_xlabel(XLabel, fontsize = FontSize_AxisLabel, fontfamily = FontName)



# Hide the right and top spines
ax.spines[['right', 'top']].set_visible(False)
ax.spines[['left', 'bottom']].set_linewidth(Axis_LineWidth)
ax.tick_params(width=Axis_LineWidth)


# save
if Save_Fig:
    
    # Create Dir
    CreateDestinationFolder(Destination_Path = Save_Dir)
    
    save_path = os.path.join(Save_Dir, Save_Name)
    
    plt.savefig(save_path, dpi='figure')    
    
    