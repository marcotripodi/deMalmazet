# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 13:53:51 2024

@author: Daniel
"""



import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


#%%
DF_SaveDir = r"O:\OneDrive - MRC Laboratory of Molecular Biology\Work\Projects\Ongoing\Inhibition\Experiments\Completed\Imaging\All\Code\motor tuned neurons"

df_All_Gyro = pd.read_pickle(os.path.join(DF_SaveDir, "df_All_Gyro"))


#%%
Thres_PValue        = .05
Thres_PearsonCoeff  = 0.1 # .1
Thres_ZScore        = 5

# CCK
gene = "CCK"

# Get # mice, recordings and neurons
nNeurons_CCK = len(df_All_Gyro.loc[(df_All_Gyro["gene"] == gene)])

Mouse_rec = [mouse_rec_neuron[:mouse_rec_neuron.rfind("_")] for mouse_rec_neuron in df_All_Gyro.loc[(df_All_Gyro["gene"] == gene)].index]

Recs = np.unique(Mouse_rec)

Mice = np.unique([mouse[:mouse.find("_")] for mouse in Recs])

print("Recorded %i %s x VGAT neurons in %i recordings in %i mice" % (nNeurons_CCK, gene, len(Recs), len(Mice)))


MotorTunedNeurons_CCK = df_All_Gyro.loc[(df_All_Gyro["gene"] == gene) * \
                            (df_All_Gyro["ZScore_Max"] >= Thres_ZScore) * \
                            (df_All_Gyro["PearsonCoeff_absGyro"].abs() >= Thres_PearsonCoeff) * \
                            (df_All_Gyro["PearsonCoeff_absGyro_PValue"].abs() <= Thres_PValue)
                            ]
     
print("Percentage motor tuned neurons in %s x VGAT: %.02f\n" % 
      (gene,
       len(MotorTunedNeurons_CCK) / sum(df_All_Gyro["gene"] == gene)))


# SST
gene = "SST"

# Get # mice, recordings and neurons
nNeurons_CCK = len(df_All_Gyro.loc[(df_All_Gyro["gene"] == gene)])

Mouse_rec = [mouse_rec_neuron[:mouse_rec_neuron.rfind("_")] for mouse_rec_neuron in df_All_Gyro.loc[(df_All_Gyro["gene"] == gene)].index]

Recs = np.unique(Mouse_rec)

Mice = np.unique([mouse[:mouse.find("_")] for mouse in Recs])

print("Recorded %i %s x VGAT neurons in %i recordings in %i mice" % (nNeurons_CCK, gene, len(Recs), len(Mice)))


MotorTunedNeurons_SST = df_All_Gyro.loc[(df_All_Gyro["gene"] == gene) * \
                            (df_All_Gyro["ZScore_Max"] >= Thres_ZScore) * \
                            (df_All_Gyro["PearsonCoeff_absGyro"].abs() >= Thres_PearsonCoeff) * \
                            (df_All_Gyro["PearsonCoeff_absGyro_PValue"].abs() <= Thres_PValue)
                            ]
                            
print("Percentage motor tuned neurons in %s x VGAT: %.02f\n" % 
      (gene,
       len(MotorTunedNeurons_SST) / sum(df_All_Gyro["gene"] == gene)))


# PV
gene = "PV"

# Get # mice, recordings and neurons
nNeurons_CCK = len(df_All_Gyro.loc[(df_All_Gyro["gene"] == gene)])

Mouse_rec = [mouse_rec_neuron[:mouse_rec_neuron.rfind("_")] for mouse_rec_neuron in df_All_Gyro.loc[(df_All_Gyro["gene"] == gene)].index]

Recs = np.unique(Mouse_rec)

Mice = np.unique([mouse[:mouse.find("_")] for mouse in Recs])

print("Recorded %i %s x VGAT neurons in %i recordings in %i mice" % (nNeurons_CCK, gene, len(Recs), len(Mice)))


MotorTunedNeurons_PV = df_All_Gyro.loc[(df_All_Gyro["gene"] == gene) * \
                            (df_All_Gyro["ZScore_Max"] >= Thres_ZScore) * \
                            (df_All_Gyro["PearsonCoeff_absGyro"].abs() >= Thres_PearsonCoeff) * \
                            (df_All_Gyro["PearsonCoeff_absGyro_PValue"].abs() <= Thres_PValue)
                            ]

               
# MotorTunedNeurons[["PearsonCoeff_absGyro_PValue", "PearsonCoeff_Gyro_PValue",
#                    "PearsonCoeff_absGyro", "PearsonCoeff_Gyro"]]


print("Percentage motor tuned neurons in %s x VGAT: %.02f\n" % 
      (gene,
       len(MotorTunedNeurons_PV) / sum(df_All_Gyro["gene"] == gene)))
    

# df_All_Gyro.loc[(df_All_Gyro["gene"] == gene) * \
#                 (df_All_Gyro["ZScore_Max"] >= Thres_ZScore)]


#%%

print("PearsonCoeff_absGyro SST: Median: %.2f, Interquartile range: [%.2f, %.2f]" % (np.percentile(MotorTunedNeurons_SST["PearsonCoeff_absGyro"], 50), 
                                                                       np.percentile(MotorTunedNeurons_SST["PearsonCoeff_absGyro"], 25),
                                                                       np.percentile(MotorTunedNeurons_SST["PearsonCoeff_absGyro"], 75)))

print("# Motor tuned neurons SST: %i\n" % (len(MotorTunedNeurons_SST["PearsonCoeff_absGyro"])))


print("PearsonCoeff_absGyro CCK: %.2f, Interquartile range: [%.2f, %.2f]" % (np.percentile(MotorTunedNeurons_CCK["PearsonCoeff_absGyro"], 50), 
                                                                          np.percentile(MotorTunedNeurons_CCK["PearsonCoeff_absGyro"], 25),
                                                                          np.percentile(MotorTunedNeurons_CCK["PearsonCoeff_absGyro"], 75)))

print("# Motor tuned neurons CCK: %i\n" % (len(MotorTunedNeurons_CCK["PearsonCoeff_absGyro"])))


print("PearsonCoeff_absGyro PV: %.2f, Interquartile range: [%.2f, %.2f]" % (np.percentile(MotorTunedNeurons_PV["PearsonCoeff_absGyro"], 50), 
                                                                          np.percentile(MotorTunedNeurons_PV["PearsonCoeff_absGyro"], 25),
                                                                          np.percentile(MotorTunedNeurons_PV["PearsonCoeff_absGyro"], 75)))

print("# Motor tuned neurons PV: %i\n" % (len(MotorTunedNeurons_PV["PearsonCoeff_absGyro"])))


#%% Box plot Coeff Pearson Abs Gyro
save_flag   = False
save_dir    = r"O:\OneDrive - MRC Laboratory of Molecular Biology\General - VGAT Paper\Figures\Paper\Figure3\Ver10"
save_name   = "BoxPlot_CorrCoeff_absgyro.svg"



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

YLim        = [-.39, .43]
YTicks      = [-0.25, 0, .25]
YLabel      = "Pearson Coeff. with Gyro"


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
                        MotorTunedNeurons_SST["PearsonCoeff_absGyro"], 
                        MotorTunedNeurons_CCK["PearsonCoeff_absGyro"], 
                        MotorTunedNeurons_PV["PearsonCoeff_absGyro"], 
                        # MotorTunedNeurons_SST["PearsonCoeff_Gyro"], 
                        # MotorTunedNeurons_CCK["PearsonCoeff_Gyro"], 
                        # MotorTunedNeurons_PV["PearsonCoeff_Gyro"], 
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
