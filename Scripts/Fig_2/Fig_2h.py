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
DF_SaveDir = r""

df_All_Gyro = pd.read_pickle(os.path.join(DF_SaveDir, "df_All_Gyro"))


#%%
Thres_3alue        = .05
Thres_PearsonCoeff  = 0.1 # .1
Thres_ZScore        = 5

# 2
gene = "2"

# Get # mice, recordings and neurons
nNeurons_2 = len(df_All_Gyro.loc[(df_All_Gyro["gene"] == gene)])

Mouse_rec = [mouse_rec_neuron[:mouse_rec_neuron.rfind("_")] for mouse_rec_neuron in df_All_Gyro.loc[(df_All_Gyro["gene"] == gene)].index]

Recs = np.unique(Mouse_rec)

Mice = np.unique([mouse[:mouse.find("_")] for mouse in Recs])

print("Recorded %i %s x VGAT neurons in %i recordings in %i mice" % (nNeurons_2, gene, len(Recs), len(Mice)))


MotorTunedNeurons_2 = df_All_Gyro.loc[(df_All_Gyro["gene"] == gene) * \
                            (df_All_Gyro["ZScore_Max"] >= Thres_ZScore) * \
                            (df_All_Gyro["PearsonCoeff_absGyro"].abs() >= Thres_PearsonCoeff) * \
                            (df_All_Gyro["PearsonCoeff_absGyro_3alue"].abs() <= Thres_3alue)
                            ]
     
print("Percentage motor tuned neurons in %s x VGAT: %.02f\n" % 
      (gene,
       len(MotorTunedNeurons_2) / sum(df_All_Gyro["gene"] == gene)))


# 1
gene = "1"

# Get # mice, recordings and neurons
nNeurons_2 = len(df_All_Gyro.loc[(df_All_Gyro["gene"] == gene)])

Mouse_rec = [mouse_rec_neuron[:mouse_rec_neuron.rfind("_")] for mouse_rec_neuron in df_All_Gyro.loc[(df_All_Gyro["gene"] == gene)].index]

Recs = np.unique(Mouse_rec)

Mice = np.unique([mouse[:mouse.find("_")] for mouse in Recs])

print("Recorded %i %s x VGAT neurons in %i recordings in %i mice" % (nNeurons_2, gene, len(Recs), len(Mice)))


MotorTunedNeurons_1 = df_All_Gyro.loc[(df_All_Gyro["gene"] == gene) * \
                            (df_All_Gyro["ZScore_Max"] >= Thres_ZScore) * \
                            (df_All_Gyro["PearsonCoeff_absGyro"].abs() >= Thres_PearsonCoeff) * \
                            (df_All_Gyro["PearsonCoeff_absGyro_3alue"].abs() <= Thres_3alue)
                            ]
                            
print("Percentage motor tuned neurons in %s x VGAT: %.02f\n" % 
      (gene,
       len(MotorTunedNeurons_1) / sum(df_All_Gyro["gene"] == gene)))


# 3
gene = "3"

# Get # mice, recordings and neurons
nNeurons_2 = len(df_All_Gyro.loc[(df_All_Gyro["gene"] == gene)])

Mouse_rec = [mouse_rec_neuron[:mouse_rec_neuron.rfind("_")] for mouse_rec_neuron in df_All_Gyro.loc[(df_All_Gyro["gene"] == gene)].index]

Recs = np.unique(Mouse_rec)

Mice = np.unique([mouse[:mouse.find("_")] for mouse in Recs])

print("Recorded %i %s x VGAT neurons in %i recordings in %i mice" % (nNeurons_2, gene, len(Recs), len(Mice)))


MotorTunedNeurons_3 = df_All_Gyro.loc[(df_All_Gyro["gene"] == gene) * \
                            (df_All_Gyro["ZScore_Max"] >= Thres_ZScore) * \
                            (df_All_Gyro["PearsonCoeff_absGyro"].abs() >= Thres_PearsonCoeff) * \
                            (df_All_Gyro["PearsonCoeff_absGyro_3alue"].abs() <= Thres_3alue)
                            ]

               
# MotorTunedNeurons[["PearsonCoeff_absGyro_3alue", "PearsonCoeff_Gyro_3alue",
#                    "PearsonCoeff_absGyro", "PearsonCoeff_Gyro"]]


print("Percentage motor tuned neurons in %s x VGAT: %.02f\n" % 
      (gene,
       len(MotorTunedNeurons_3) / sum(df_All_Gyro["gene"] == gene)))
    

# df_All_Gyro.loc[(df_All_Gyro["gene"] == gene) * \
#                 (df_All_Gyro["ZScore_Max"] >= Thres_ZScore)]


#%%

print("PearsonCoeff_absGyro 1: Median: %.2f, Interquartile range: [%.2f, %.2f]" % (np.percentile(MotorTunedNeurons_1["PearsonCoeff_absGyro"], 50), 
                                                                       np.percentile(MotorTunedNeurons_1["PearsonCoeff_absGyro"], 25),
                                                                       np.percentile(MotorTunedNeurons_1["PearsonCoeff_absGyro"], 75)))

print("# Motor tuned neurons 1: %i\n" % (len(MotorTunedNeurons_1["PearsonCoeff_absGyro"])))


print("PearsonCoeff_absGyro 2: %.2f, Interquartile range: [%.2f, %.2f]" % (np.percentile(MotorTunedNeurons_2["PearsonCoeff_absGyro"], 50), 
                                                                          np.percentile(MotorTunedNeurons_2["PearsonCoeff_absGyro"], 25),
                                                                          np.percentile(MotorTunedNeurons_2["PearsonCoeff_absGyro"], 75)))

print("# Motor tuned neurons 2: %i\n" % (len(MotorTunedNeurons_2["PearsonCoeff_absGyro"])))


print("PearsonCoeff_absGyro 3: %.2f, Interquartile range: [%.2f, %.2f]" % (np.percentile(MotorTunedNeurons_3["PearsonCoeff_absGyro"], 50), 
                                                                          np.percentile(MotorTunedNeurons_3["PearsonCoeff_absGyro"], 25),
                                                                          np.percentile(MotorTunedNeurons_3["PearsonCoeff_absGyro"], 75)))

print("# Motor tuned neurons 3: %i\n" % (len(MotorTunedNeurons_3["PearsonCoeff_absGyro"])))


#%% Box plot Coeff Pearson Abs Gyro
save_flag   = False
save_dir    = r""
save_name   = "BoxPlot_CorrCoeff_absgyro.svg"



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
                        MotorTunedNeurons_1["PearsonCoeff_absGyro"], 
                        MotorTunedNeurons_2["PearsonCoeff_absGyro"], 
                        MotorTunedNeurons_3["PearsonCoeff_absGyro"], 
                        # MotorTunedNeurons_1["PearsonCoeff_Gyro"], 
                        # MotorTunedNeurons_2["PearsonCoeff_Gyro"], 
                        # MotorTunedNeurons_3["PearsonCoeff_Gyro"], 
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
