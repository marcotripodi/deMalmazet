# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 14:49:45 2026

@author: Daniel
"""


import os
import numpy as np
import pandas as pd
from cellpose import io
from tqdm import tqdm
import matplotlib.pyplot as plt
from PathGeneral_Func import CreateDestinationFolder


#%% Set up paths and data files
Data_Dir = r"O:\OneDrive - MRC Laboratory of Molecular Biology\Work\Projects\Ongoing\Inhibition\Experiments\Completed\GQ Dreadd validation with cfos\Data"

SegFiles_Dir = os.path.join(Data_Dir, 
                            "Seg_Files")

CFOS_Images_Dir = os.path.join(Data_Dir, 
                               "CFOS") 

ROIs_Dir = os.path.join(Data_Dir, 
                        "ROIs") 


Exp_Names = ["MTBZ11_1a_Slide2_slice09",
             "MTBZ11_1a_Slide3_slice06",
             "MTBZ11_1a_Slide4_slice07",
            
             "MTBZ11_1b_Slide3_slice04",
             "MTBZ11_1b_Slide3_slice14",
             "MTBZ11_1b_Slide3_slice04",
            
             "MTBZ11_1c_Slide2_slice10",
             "MTBZ11_1c_Slide3_slice05",
             "MTBZ11_1c_Slide3_slice13",
                          
             "MTBZ11_1d_Slide3_slice09",
             "MTBZ11_1d_Slide3_slice11",
             "MTBZ11_1d_Slide3_slice14",
            ]


#%% Load
df_Main = pd.DataFrame()
for exp_name in tqdm(Exp_Names):
# for exp_name in [Exp_Names[0]]:
    
    print("\n\n\nProcessing exp %s" % (exp_name))

    Seg_Filename = r"mCherry_GammaCorrected_" + exp_name + "_seg.npy"
    CFOS_Filename = exp_name + r".tif"
    
    
    # Load image and mask
    dat = np.load(os.path.join(SegFiles_Dir, Seg_Filename), allow_pickle=True).item()
    img = io.imread(os.path.join(CFOS_Images_Dir, CFOS_Filename))


    # Get number of masks
    nROIs = np.max(np.unique(dat['masks']))
    print("\n%i ROIs detected\n\n" % (nROIs))
    
    
    #% Calculate max and mean pixel value for each roi
    ROIs_Mean   = np.empty(nROIs)
    ROIs_Max    = np.empty(nROIs)    
    for roi in range(1, nROIs + 1):
        
        ROI_Pixels = img[dat['masks'] == roi]
        
        # ROIs_Mean[roi-1] = np.mean(ROI_Pixels)
        ROIs_Max[roi-1] = np.max(ROI_Pixels)
        
    
    # Compute threshold for considering cell CFOS positive based on the background value
    Threshold_CFOSpos = np.mean(img[dat['masks'] == 0]) + np.std(img[dat['masks'] == 0]) * 3
    # Threshold_CFOSpos = np.mean(img) + np.std(img) * 3
    
    
    # Get number CFOS positive
    nCells_pos = np.sum(ROIs_Max > Threshold_CFOSpos)
    nCells_neg = np.sum(ROIs_Max <= Threshold_CFOSpos)
    
    
    # Create dataframe
    if exp_name[8] == "a" or exp_name[8] == "b":
        Injection = "DCZ"
    else:
        Injection = "PBS"
        
        
    df = pd.DataFrame(data = {"Injection" : Injection, 
                              "Total cells" : nROIs,
                              "CFOS pos cells" : nCells_pos,
                              "% positive" : nCells_pos/nROIs * 100}, 
                      index = [exp_name])
    
    df_Main = pd.concat([df_Main, df])
    
    
#%% plot    
save_fig    = False
save_dir    = r"O:\OneDrive - MRC Laboratory of Molecular Biology\General - VGAT Paper\Figures\Paper\Figure6\Ver6"
save_name   = "BarPlot_CFOSpos.svg"


# Param figure
BarWidths               = .5
Space_Bars_withinGene   = .4
Space_Bars_accrossGene  = .8
Xticklable_angle        = 15



FontSize_TickLabel  = 6
FontSize_AxisLabel  = 7
FontName            = "arial"


Color = "k"
Color = ".3"

Ticks = np.cumsum([Space_Bars_withinGene, Space_Bars_accrossGene])
    
    
# Error Bars
capsize     = 6
capthick    = .6
elinewidth  = .6


XTicks      = Ticks
XTickLabels = ["DCZ", "PBS"]

XLabel = ""

YTicks  = [0, 10, 20, 30]
# YLabel  = "% of c-Fos positive cells \n among mCherry positive cells"
YLabel  = ""
YLim    = [0, 40.5]
# YLim    = [0, 3.5+df_Main["% positive"][df_Main["Injection"] == "DCZ"].mean() + df_Main["% positive"][df_Main["Injection"] == "PBS"].std()/2]

# figure
Fig_Width = 1.2
Fig_Height = 1


# Axis
Axis_left       = .3
Axis_bottom     = 0.2
Axis_width      = .6
Axis_height     = .7
Axis_LineWidth  = .5



#% Plot the PDF.
Fig = plt.figure(figsize = [Fig_Width, Fig_Height], dpi = 300)

ax = plt.axes([Axis_left, Axis_bottom, Axis_width, Axis_height])   


plt.bar(Ticks[0], 
        df_Main["% positive"][df_Main["Injection"] == "DCZ"].mean(),
        color = Color,
        width = BarWidths,
        yerr = df_Main["% positive"][df_Main["Injection"] == "DCZ"].std()/2,
        ecolor = Color,
        capsize = capsize,
        error_kw = {"elinewidth": elinewidth,
                    "capthick": capthick})


plt.bar(Ticks[1], 
        df_Main["% positive"][df_Main["Injection"] == "PBS"].mean(),
        color = Color,
        width = BarWidths,
        yerr = df_Main["% positive"][df_Main["Injection"] == "PBS"].std()/2,
        ecolor = Color,
        capsize = capsize,
        error_kw = {"elinewidth": elinewidth,
                    "capthick": capthick})



ax.set_xlabel(XLabel, fontsize = FontSize_AxisLabel, fontfamily = FontName)
ax.set_xlim([Ticks[0] - Space_Bars_withinGene, Ticks[-1] + Space_Bars_withinGene])
ax.set_xticks(XTicks)
ax.set_xticklabels(XTickLabels, 
                   fontsize = FontSize_TickLabel, 
                   fontfamily = FontName)


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
        
    
#%% Stats
print("\nDCZ: mean +- std:\n%.2f +- %.2f" % (df_Main["% positive"][df_Main["Injection"] == "DCZ"].mean(), 
                                            df_Main["% positive"][df_Main["Injection"] == "DCZ"].std()))
    
print("\nPBS: mean +- std:\n%.2f +- %.2f" % (df_Main["% positive"][df_Main["Injection"] == "PBS"].mean(), 
                                            df_Main["% positive"][df_Main["Injection"] == "PBS"].std()))
    