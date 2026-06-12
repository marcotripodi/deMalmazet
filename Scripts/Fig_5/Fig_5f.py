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
import seaborn as sns



#%%
nSlices = 5

df_Main = pd.DataFrame()

for Slice_Num in range(1, nSlices + 1):

    print("\n\n\nProcessing slice %s" % (Slice_Num))
    
    # Set up path
    Seg_Dir     = r"\Slice_" + str(Slice_Num) + "\Corrected"
    Image_Dir   = r"\Slice_" + str(Slice_Num) + "\Cropped"

    
    #% Load
    Exp_Name = "MTBZ9.4b_Sd4Sc" + str(Slice_Num) + "_"
    
    
    Seg_Filename    = Exp_Name + r"GFP_FarRed_seg.npy"
    VGAT_Filename   = Exp_Name + r"VGAT_GFP.tif"
    Marker2_Filename    = Exp_Name + r"2_Cy3.tif"
    
    
    dat         = np.load(os.path.join(Seg_Dir, Seg_Filename), allow_pickle=True).item()
    img_VGAT    = io.imread(os.path.join(Image_Dir, VGAT_Filename))
    img_2     = io.imread(os.path.join(Image_Dir, Marker2_Filename))

    
    #% Calculate max and mean pixel value for each roi
    nGFP_Pos_Cells = np.max(np.unique(dat['masks']))
    print("\n%i ROIs detected" % (nGFP_Pos_Cells))
    
    ROIs_Max_2    = np.empty(nGFP_Pos_Cells)
    ROIs_Max_VGAT   = np.empty(nGFP_Pos_Cells)
    for roi in tqdm(range(1, nGFP_Pos_Cells + 1)):
    
        # Get Roi pixels
        ROI_Pixels = img_VGAT[dat['masks'] == roi]
        
        ROIs_Max_VGAT[roi-1] = np.max(ROI_Pixels)
    
    
        # Get Roi pixels
        ROI_Pixels = img_2[dat['masks'] == roi]
        
        ROIs_Max_2[roi-1] = np.max(ROI_Pixels)
    
        
    #     
    # Threshold_2   = 3000
    # Threshold_2   = np.mean(img_2) + 3*np.std(img_2)
    Threshold_2   = np.mean(img_2[dat['masks'] == 0]) + 3*np.std(img_2[dat['masks'] == 0])
    # Threshold_2   = np.mean(img_2[dat['masks'] == 0]) + 4*np.std(img_2[dat['masks'] == 0])
    
    # Threshold_VGAT  = 6000
    # Threshold_VGAT  = np.mean(img_VGAT) + 3*np.std(img_VGAT)
    Threshold_VGAT  = np.mean(img_VGAT[dat['masks'] == 0]) + 3*np.std(img_VGAT[dat['masks'] == 0])
    # Threshold_VGAT  = np.mean(img_VGAT[dat['masks'] == 0]) + 4*np.std(img_VGAT[dat['masks'] == 0])
    
    Marker2_VGAT_Cells = np.flatnonzero((ROIs_Max_VGAT > Threshold_VGAT) * (ROIs_Max_2 > Threshold_2))
    
    
    print("\nPercentage 2 VGAT neurons among GFP pos cell:", len(Marker2_VGAT_Cells) / nGFP_Pos_Cells * 100)

    
        
    df = pd.DataFrame(data = {
                              "Total gfp cells" : nGFP_Pos_Cells,
                              "# 2_VGAT_Cells" : len(Marker2_VGAT_Cells),
                              "% 2_VGAT_Cells" : len(Marker2_VGAT_Cells) / nGFP_Pos_Cells * 100}, 
                      index = [Slice_Num])
    
    df_Main = pd.concat([df_Main, df])
    
    

#%% Strip_Plot
save_fig    = False
save_dir    = r""
save_name   = "StriPlot_percentage_2_VGAT_Cells.svg"


MarkerSize = 2.3


# palette = sns.color_palette("colorblind", n_colors = 2)
palette = sns.color_palette("colorblind", n_colors = 1)

Color_1   = [141/255, 160/255, 203/255]
Color_3    = [252/255, 141/255, 98/255]
Color_2   = [102/255, 194/255, 165/255]


# Line_PVue
Line_PVue_Color       = [.6, .6, .6]
Line_PVue_LineWidth   = .6
Line_PVue_LineStyle   = "--"


# X axis
XTicks      = [0]
xticklabels = ["$\mathregular{2^{ON}xVGAT^{ON}}$"]
xticklabels = [""]
XLabel      = ""
XLim        = [-.1, .1]


# Y axis
yslack          = 2
YLabel          = "% of $\mathregular{2^{ON}xVGAT^{ON}}$ \n among gfp positive neurons"
yticks          = [0, 30, 60]
yticks_label    = yticks #["0", "10"]
YLim            = [0, df_Main["% 2_VGAT_Cells"].max() + yslack]


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

ErrorBar_LineWidth = .8
ErrorBar_Color = "k"
ErrorBar_XPos = -0.05

Fig = plt.figure(figsize = [Fig_Width, Fig_Height], dpi = 300)

ax = plt.axes([Axis_left, Axis_bottom, Axis_width, Axis_height])   
    

ax = sns.swarmplot( ax = ax, 
                    data = df_Main,
                    
                    y = "% 2_VGAT_Cells",
                    hue = None,
                    color = ".3",
                    size = MarkerSize,
                    linewidth = Axis_LineWidth,
                    legend = False, 
                    edgecolor = "none")


# Mean error bars
ErrorBar = ax.errorbar(ErrorBar_XPos,
                        df_Main["% 2_VGAT_Cells"].mean(),
                        marker = "",
                        c = ErrorBar_Color,
                        linewidth = ErrorBar_LineWidth,
                        yerr = df_Main["% 2_VGAT_Cells"].std() / 2,
                        ecolor = ErrorBar_Color,
                        elinewidth = ErrorBar_LineWidth,
                        capsize = 2.5,
                        capthick = ErrorBar_LineWidth,)



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
if save_fig:
    
    # Create Dir
    CreateDestinationFolder(Destination_Path = save_dir)
    
    save_path = os.path.join(save_dir, save_name)
    
    plt.savefig(save_path, dpi='figure')    
    
    
#%%
print("\n\nMean +- STD:")
print("%.2f +- %.2f" % (df_Main["% 2_VGAT_Cells"].mean(), 
                        df_Main["% 2_VGAT_Cells"].std())
      )


print("\n\nMean +- STD:")
print("%.2f +- %.2f" % (df_Main["% 2_VGAT_Cells"][:-1].mean(), 
                        df_Main["% 2_VGAT_Cells"][:-1].std())
      )

            