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


#%%
DB_dir  = r"O:\OneDrive - MRC Laboratory of Molecular Biology\Work\Projects\Ongoing\Inhibition\Experiments\Completed\Imaging\SST\Data\DB"

DB_name = "DB_VGAT_SST_SingleRecs_Dark_Smooth"


with open(os.path.join(DB_dir, DB_name), 'rb') as handle:
    DB_Recs = pickle.load(handle)
    
list(DB_Recs['MTBZ5.3a_20230906'])    
#%%
RecsOfInterest = 'MTBZ5.3a_20230906'


# Set Save dir
SaveDir = CreateDestinationFolder(
    Destination_Path = r"O:\OneDrive - MRC Laboratory of Molecular Biology\General - VGAT Paper\Figures\Paper\Figure3\Ver3\Example_Traces")

SaveName    = RecsOfInterest + "_withHistogram.svg"
SaveFig     = False

Frames2Plot = [2000, 12000]

# Figure Params
Fig_Width   = 2.6
Fig_Height  = 1.2

Color_SST   = "#7570b3"
Color_Gyro  = (1, 0.59, 0)

Color_Traces = [Color_SST]
Trace_LineWidth = .7


# Axis
Axis_left       = .16
Axis_bottom     = 0.05
Axis_right      = .79
Axis_top        = .94
Axis_LineWidth  = .5


# Y Axis
YTicks_DFoverF  = [-.5, 0, .5]
YLabel_DFoverF  = "DF over F"

YTicks_Gyro     = [-100, 100]
YLabel_Gyro     = "gyro (\u00b0/s)"

# X axis
XTicks          = []
Offset_YAxis    = 5


# Font params
FontSize_TickLabel = 6
FontSize_AxisLabel = 7
FontName = "arial"


# Heat map params
Ticks_LineWidth = .2
Ticks_Length    = 3


#% Histograms axis params
Axis_Hist_offset    = .05
Axis_Hist_width     = .15


# Param figure
BinSize         = .05
rwidth          = 1
alpha           = .6
density_flag    = True
Color_Hist      = Color_SST


# X axis
XTicks_Hist = [0, 1]
XLabel_hist = "Probability\ndensity"


#% middle line
MiddleLine_Color = ".5" # [.6, .6, 6]
MiddleLine_LineWidth = .5
MiddleLine_LineStyle = "--"


# XScale
XScale_Pos0         = 1000
XScale_time         = 20 # in seconds
XScale_Linewidth    = 2
XScale_Color        = "0"
XScale_Text         = str(XScale_time) + "s"
XScale_Text_Y       = 25
XScale_Text_X       = 0


fig, ax = plt.subplots(2, 1, 
                       figsize = (Fig_Width, Fig_Height),
                       sharex = True, 
                       dpi = 300,
                       gridspec_kw = {"left": Axis_left, 
                                      "right": Axis_right, 
                                      "top": Axis_top, 
                                      "bottom": Axis_bottom})

              
# Plot DF trace
trace = DB_Recs[RecsOfInterest]["DFoverF"][0, Frames2Plot[0]:Frames2Plot[1]]

ax[0].plot(trace,
                c = Color_Traces[0],
                linewidth = Trace_LineWidth)


# Yaxis
YLim = [np.min(trace), np.max(trace)]
ax[0].set_ylim(YLim)
ax[0].set_yticks(YTicks_DFoverF)
ax[0].set_yticklabels(YTicks_DFoverF, fontsize = FontSize_TickLabel, fontfamily = FontName)

# Ylabel
ax[0].set_ylabel(YLabel_DFoverF, fontsize = FontSize_AxisLabel, fontfamily = FontName)



# Plot hist        
        
# set Axis his
Bbox = ax[0].get_position()
Axis_Hist_left      = Bbox.get_points()[1, 0] + Axis_Hist_offset        
Axis_Hist_bottom    = Bbox.get_points()[0, 1]
Axis_Hist_height    = Bbox.bounds[-1]

ax_Hist = plt.axes([Axis_Hist_left, Axis_Hist_bottom, Axis_Hist_width, Axis_Hist_height])   


Bin_Edges = np.arange(np.min(trace), np.max(trace) + BinSize, BinSize)

Bin_counts, _, _ = ax_Hist.hist(trace, 
                                 bins = Bin_Edges, 
                                 density = density_flag, 
                                 alpha=alpha, 
                                 color = Color_Traces[0], 
                                 rwidth = rwidth, 
                                 orientation=u'horizontal')


# y axis
ax_Hist.set_ylim(YLim)
ax_Hist.set_yticks(YTicks_DFoverF)
# ax_Hist.set_yticklabels(YTicks_DFoverF, fontsize = FontSize_TickLabel, fontfamily = FontName)
ax_Hist.set_yticklabels(["" for i in YTicks_DFoverF])


# X axis
Xlim_hist = [0, Bin_counts.max()]
ax_Hist.set_xlim(Xlim_hist)
ax_Hist.set_xticks(XTicks_Hist)
ax_Hist.set_xticklabels(["" for i in XTicks_Hist])


ax_Hist.set_xticklabels(XTicks_Hist, fontsize = FontSize_TickLabel, fontfamily = FontName)
ax_Hist.set_xlabel(XLabel_hist, fontsize = FontSize_AxisLabel, fontfamily = FontName)


# Plot middle line
Middle = (YLim[1] + YLim[0])/2
ax_Hist.plot([0, Xlim_hist[1]], 
             [Middle, Middle],
             c = MiddleLine_Color,
             linewidth = MiddleLine_LineWidth,
             linestyle = MiddleLine_LineStyle)
 
    
# Remove spines
ax_Hist.spines[['top', "right"]].set_visible(False)
ax_Hist.spines[['left', 'bottom']].set_linewidth(Axis_LineWidth)

# increase tick width
ax_Hist.tick_params(size = Ticks_Length, width = Ticks_LineWidth)

    

# Gyro
ax[1].plot(DB_Recs[RecsOfInterest]["Gyro_Smooth"][Frames2Plot[0]:Frames2Plot[1]],
           c = Color_Gyro,
           linewidth = Trace_LineWidth)

# Y axis
YLim_Gyro = [np.min(DB_Recs[RecsOfInterest]["Gyro_Smooth"][Frames2Plot[0]:Frames2Plot[1]]),
             np.max(DB_Recs[RecsOfInterest]["Gyro_Smooth"][Frames2Plot[0]:Frames2Plot[1]])
             ]
ax[1].set_ylim(YLim_Gyro)
ax[1].set_yticks(YTicks_Gyro)
ax[1].set_yticklabels(YTicks_Gyro, fontsize = FontSize_TickLabel, fontfamily = FontName)
ax[1].set_ylabel(YLabel_Gyro, fontsize = FontSize_AxisLabel, fontfamily = FontName)

            
for ax_subplot in ax:
    ax_subplot.spines[['right', 'top', 'bottom']].set_visible(False)
    ax_subplot.spines['left'].set_linewidth(Axis_LineWidth)
    ax_subplot.tick_params(size = Ticks_Length, width = Ticks_LineWidth)
    # ax_subplot.axis("off")
    

# X axis
ax[1].set_xticks(XTicks)
# ax.set_xlabel("Timepoints", fontsize = FontSize_AxisLabel, fontfamily = FontName)
# ax.set_xticklabels(ax.get_xticks(), fontsize = FontSize_TickLabel, fontfamily = FontName)


nFrames2Plot = len(DB_Recs[RecsOfInterest]["DFoverF"][0, Frames2Plot[0]:Frames2Plot[1]])
ax[1].set_xlim([-Offset_YAxis, nFrames2Plot])
            
    
# Add time scale
XScale_Frames = XScale_time * DB_Recs[RecsOfInterest]["FrameRate"]
ax[1].plot([XScale_Pos0, XScale_Pos0 + XScale_Frames],
           [YLim_Gyro[0], YLim_Gyro[0]],
           linewidth = XScale_Linewidth,
           c = XScale_Color)
           
# Add text
ax[1].text(XScale_Pos0 + XScale_Text_X,
           YLim_Gyro[0] + XScale_Text_Y,
           s = XScale_Text,
           fontsize = FontSize_TickLabel, 
           fontfamily = FontName
           )

# save
if SaveFig:
    plt.savefig(os.path.join(SaveDir, SaveName), dpi='figure')
        