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


# plt.close("all")
#%%
DB_dir  = r""

DB_name = ""


with open(os.path.join(DB_dir, DB_name), 'rb') as handle:
    DB_Recs = pickle.load(handle)
    
    

RecsOfInterest = 'MTCE1.3a_20210518_z53'

Neuron_Tonic = 3

# for Neuron_Phasic in [8,25,72,76,88,89]: #25, 89
for Neuron_Phasic in [25]: # [8,11,25,26,31,32,52,55,56,58,72,76,88,89,100, 101]


    Neurons2Plot = [Neuron_Tonic, Neuron_Phasic]
        
    # Set Save dir
    SaveDir = CreateDestinationFolder(
        Destination_Path = r"")
    
    SaveName    = RecsOfInterest + "_" + str(Neuron_Tonic) + "_" + str(Neuron_Phasic) + "_withHistogram.svg"
    SaveFig     = False
    
    Frames2Plot = [0, 30000]
    
    # Figure Params
    Fig_Width   = 2.6
    Fig_Height  = 1.2
    
    Color_Tonic     = (0, 0.68, 0.94)
    Color_Phasic    = (0.93, 0, 0.55)
    Color_Gyro      = (1, 0.59, 0)
    
    Color_Traces = [Color_Tonic, Color_Phasic]
    Trace_LineWidth = .7
    
    
    # Axis
    Axis_left       = .14
    Axis_bottom     = 0.05
    Axis_right      = .79
    Axis_top        = .94
    Axis_LineWidth  = .5
    
    
    # Y Axis
    YTicks_DFoverF  = [0, 1]
    YLabel_DFoverF  = "DF over F"
    
    YTicks_Gyro     = [-50, 0]
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
    Color_Hist      = Color_Tonic
    
    
    # X axis
    XTicks_Hist = [0, 1]
    XLabel_hist = "Probability\ndensity"
    
    
    #% middle line
    Thres_Ratio             = 1/2
    MiddleLine_Color = ".5" # [.6, .6, 6]
    MiddleLine_LineWidth = .5
    MiddleLine_LineStyle = "--"
    
    
    # XScale
    XScale_Pos0         = 1000
    XScale_time         = 20 # in seconds
    XScale_Linewidth    = 2
    XScale_Color        = "0"
    XScale_Text         = str(XScale_time) + "s"
    XScale_Text_Y       = 5
    XScale_Text_X       = -100
    

    fig, ax = plt.subplots(3, 1, 
                           figsize = (Fig_Width, Fig_Height),
                           sharex = True, 
                           dpi = 300,
                           label = Neuron_Phasic,
                           gridspec_kw = {"left": Axis_left, 
                                          "right": Axis_right, 
                                          "top": Axis_top, 
                                          "bottom": Axis_bottom})
    
                                                    
    #% plot traces
    df_num = 0
    for neuron2plot in Neurons2Plot:
        
        # Plot DF trace
        trace = DB_Recs[RecsOfInterest]["DFoverF"][neuron2plot, Frames2Plot[0]:Frames2Plot[1]]
        
        ax[df_num].plot(trace,
                        c = Color_Traces[df_num],
                        linewidth = Trace_LineWidth)
        
        
        # Yaxis
        YLim = [np.min(trace), np.max(trace)]
        ax[df_num].set_ylim(YLim)
        ax[df_num].set_yticks(YTicks_DFoverF)
        ax[df_num].set_yticklabels(YTicks_DFoverF, fontsize = FontSize_TickLabel, fontfamily = FontName)
        
        # Ylabel
        if df_num == 0:
            ax[df_num].set_ylabel(YLabel_DFoverF, fontsize = FontSize_AxisLabel, fontfamily = FontName)
    
    
    
        # Plot hist        
                
        # set Axis his
        Bbox = ax[df_num].get_position()
        Axis_Hist_left      = Bbox.get_points()[1, 0] + Axis_Hist_offset        
        Axis_Hist_bottom    = Bbox.get_points()[0, 1]
        Axis_Hist_height    = Bbox.bounds[-1]
        
        ax_Hist = plt.axes([Axis_Hist_left, Axis_Hist_bottom, Axis_Hist_width, Axis_Hist_height])   
        
        
        Bin_Edges = np.arange(np.min(trace), np.max(trace) + BinSize, BinSize)
        
        Bin_counts, _, _ = ax_Hist.hist(trace, 
                                         bins = Bin_Edges, 
                                         density = density_flag, 
                                         alpha=alpha, 
                                         color = Color_Traces[df_num], 
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
        
        if df_num == 1:
            ax_Hist.set_xticklabels(XTicks_Hist, fontsize = FontSize_TickLabel, fontfamily = FontName)
            ax_Hist.set_xlabel(XLabel_hist, fontsize = FontSize_AxisLabel, fontfamily = FontName)
    
    
        # Plot middle line
        Thres = (YLim[1] + YLim[0]) * Thres_Ratio
        ax_Hist.plot([0, Xlim_hist[1]], 
                     [Thres, Thres],
                     c = MiddleLine_Color,
                     linewidth = MiddleLine_LineWidth,
                     linestyle = MiddleLine_LineStyle)
     
            
        # Remove spines
        ax_Hist.spines[['top', "right"]].set_visible(False)
        ax_Hist.spines[['left', 'bottom']].set_linewidth(Axis_LineWidth)
        
        # increase tick width
        ax_Hist.tick_params(size = Ticks_Length, width = Ticks_LineWidth)
        
        
        df_num += 1
        
    
    # Gyro
    ax[2].plot(DB_Recs[RecsOfInterest]["Gyro_Smooth"][Frames2Plot[0]:Frames2Plot[1]],
               c = Color_Gyro,
               linewidth = Trace_LineWidth)
    
    # Y axis
    YLim_Gyro = [np.min(DB_Recs[RecsOfInterest]["Gyro_Smooth"][Frames2Plot[0]:Frames2Plot[1]]),
                 np.max(DB_Recs[RecsOfInterest]["Gyro_Smooth"][Frames2Plot[0]:Frames2Plot[1]])
                 ]
    ax[2].set_ylim(YLim_Gyro)
    ax[2].set_yticks(YTicks_Gyro)
    ax[2].set_yticklabels(YTicks_Gyro, fontsize = FontSize_TickLabel, fontfamily = FontName)
    ax[2].set_ylabel(YLabel_Gyro, fontsize = FontSize_AxisLabel, fontfamily = FontName)
    
                
    for ax_subplot in ax:
        ax_subplot.spines[['right', 'top', 'bottom']].set_visible(False)
        ax_subplot.spines['left'].set_linewidth(Axis_LineWidth)
        ax_subplot.tick_params(size = Ticks_Length, width = Ticks_LineWidth)
        # ax_subplot.axis("off")
        
    
    # X axis
    ax[2].set_xticks(XTicks)
    # ax.set_xlabel("Timepoints", fontsize = FontSize_AxisLabel, fontfamily = FontName)
    # ax.set_xticklabels(ax.get_xticks(), fontsize = FontSize_TickLabel, fontfamily = FontName)
    
    
    nFrames2Plot = len(DB_Recs[RecsOfInterest]["DFoverF"][neuron2plot, Frames2Plot[0]:Frames2Plot[1]])
    ax[2].set_xlim([-Offset_YAxis, nFrames2Plot])
                
        
    # Add time scale
    XScale_Frames = XScale_time * DB_Recs[RecsOfInterest]["FrameRate"]
    ax[2].plot([XScale_Pos0, XScale_Pos0 + XScale_Frames],
               [YLim_Gyro[0], YLim_Gyro[0]],
               linewidth = XScale_Linewidth,
               c = XScale_Color)
               
    # Add text
    ax[2].text(XScale_Pos0 + XScale_Text_X,
               YLim_Gyro[0] + XScale_Text_Y,
               s = XScale_Text,
               fontsize = FontSize_TickLabel, 
               fontfamily = FontName
               )

    # save
    if SaveFig:
        plt.savefig(os.path.join(SaveDir, SaveName), dpi='figure')
        