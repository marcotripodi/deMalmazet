# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 11:23:43 2024

@author: Daniel
"""



import os
import seaborn as sns
import matplotlib.pyplot as plt
from PathGeneral_Func import CreateDestinationFolder
import numpy as np


#%% Data
# cut off < 10 days

Time2Capture_Early = [
59.80718876, #MTCM2.1b
61.02327458,
103.5122721,

222.2967487, # MTBZ3.2c
147.6284136,

266.3338124, # MTCG2.2d
345.7059342,
17.62614419,
27.87654134,

]
    



Time2Capture_Late = [
17.43980795, #MTAO68.1c

35.55488312, # MTAH9.3b
38.80315501,
66.61895778,

72.57230754, #MTCM2.1b
68.13271605,
24.86584898,
28.04767712,
26.72481789,

40.73678838, #MTCM2.1a
25.0609983,
46.73301601,
90.03277055,
80.97266449,
37.72096692,
45.57833422,
45.35138815,


48.47494553, #MTCG2.2d
38.01742919,
5.938494168,

]

print("Before 10 days: mean: %.2fs, std: %.2f" % (np.mean(Time2Capture_Early), np.std(Time2Capture_Early)))
print("Adter 10 days: mean: %.2fs, std: %.2f" % (np.mean(Time2Capture_Late), np.std(Time2Capture_Late)))


#%%

Save_Dir = r"O:\OneDrive - MRC Laboratory of Molecular Biology\Work\Projects\Ongoing\Inhibition\Experiments\Completed\Cricket_HeadFixed\Analysis"
Save_Fig = False

Save_Name = "Early_Late_10days.svg"


size = 1
color = "k"
edgecolor = "k"

linewidth = 1


# X axis
XTicks      = [0, 1]
xticklabels = ["< 10", ">= 10"]
XLabel      = "Training days"
XLim        = [-.5, 1.5]


# Y axis
YLabel          = "Time to capture (s)"
yticks          = [0, 50, 100, 150, 200, 250, 300]
yticks_label    = yticks
YLim            = [0, 355]


# figure
Fig_Width   = 1.8
Fig_Height  = 1.3

    
# Axis
Axis_LineWidth  = .3
Axis_left       = .23
Axis_bottom     = 0.25
Axis_width      = .76
Axis_height     = .7


# Traces    
Linewidth_Traces    = .8
 

FontSize_TickLabel  = 6
FontSize_AxisLabel  = 7
FontName            = "arial"



Fig = plt.figure(figsize = [Fig_Width, Fig_Height], dpi = 300)

ax = plt.axes([Axis_left, Axis_bottom, Axis_width, Axis_height])   
    
sns.swarmplot(data = [Time2Capture_Early, Time2Capture_Late],
              color = color,
              edgecolor = edgecolor,
              size = size, 
              linewidth = linewidth,
              ax = ax
              )


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
ax.tick_params(width = Axis_LineWidth)


# save
if Save_Fig:
    
    # Create Dir
    CreateDestinationFolder(Destination_Path = Save_Dir)
    
    save_path = os.path.join(Save_Dir, Save_Name)
    
    plt.savefig(save_path, dpi='figure')    
    
    