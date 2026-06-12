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
import matplotlib
from sklearn.decomposition import NMF
from Prediction_Func import PlotModel3d

# plt.close("all")
#%%
DB_dir  = r"O:\OneDrive - MRC Laboratory of Molecular Biology\Work\Projects\Ongoing\Inhibition\Experiments\Completed\Imaging\VGAT\Data\DB"

DB_name = "DB_VGAT_Dark_Smooth_2"
# DB_name = "DB_VGAT_SingleRecs_Dark_Smooth"


with open(os.path.join(DB_dir, DB_name), 'rb') as handle:
    DB_Recs = pickle.load(handle)
    
 
rec = 'MTCE1.3a_20210526'


#%%
Thres_ZScore    = 5
Thres_Ratio     = .5

# NMF params
nCmpnts         = 2 #None
init            = "nndsvd" # 'random', nndsvd (better for sparseness), nndsvda (better when sparsity is not desired), nndsvdar (generally faster, less accurate alternative to NNDSVDa for when sparsity is not desired)
max_iter        = 10000
tol             = 1e-4
ErrorFunction   = "frobenius" # kullback-leibler, frobenius
solver          = "cd" # mu or cd
l1_ratio        = 1
alpha_W         = 0
alpha_H         = 0 #0.15 #'same'


# Create NMF model
model = NMF(n_components = nCmpnts, 
            init = init, 
            max_iter = max_iter,
            beta_loss = ErrorFunction,
            solver = solver,
            alpha_W = alpha_W, 
            alpha_H = alpha_H,
            l1_ratio = l1_ratio)



# Select responding neurons
Responding_Neurons = np.flatnonzero(DB_Recs[rec]["ZScore_Max"] >= Thres_ZScore)

# Get DF traces
X = DB_Recs[rec]["DFoverF"][Responding_Neurons, :].T    

nTimePoints, nNeurons = X.shape



# Divide neurons into phasic and tonic
# Get middle of traces range
Middle = (np.max(X, axis = 0) + np.min(X, axis = 0) ) / 2

# Divide timepoints into active and inactive
TimePoints_Neurons_Active   = X > Middle
TimePoints_Neurons_Inactive = X < Middle

# Compute ratio time spent active on inactive
Ratio_Active_Inactive_rec = np.sum(TimePoints_Neurons_Active, axis = 0) / np.sum(TimePoints_Neurons_Inactive, axis = 0)

# Classify neurons into tonic or phasic
Neurons_Tonic    = np.flatnonzero(Ratio_Active_Inactive_rec > Thres_Ratio)
Neurons_Phasic   = np.flatnonzero(Ratio_Active_Inactive_rec < Thres_Ratio)



# Classify neurons into go and no go
# Remove negative values from X
X_NonNeg = X.copy()
X_NonNeg[X_NonNeg < 0] = 0

#% Fit NMF
W = model.fit_transform(X_NonNeg)

# Sort components so no go is 0    
#% Compute pearson correlation
Pearson_Gyro = np.empty(nCmpnts)
for cmpnt in range(nCmpnts):
    
    # Correlation with absolute gyro
    Pearson_Gyro[cmpnt] = np.corrcoef(W[:, cmpnt], np.abs(DB_Recs[rec]["Gyro_Smooth"]))[0, 1]


# Sort so No go componenet is first
H = model.components_[np.argsort(Pearson_Gyro), :]


# Get phasic go and no go neurons
Neurons_PhasicNoGO  = np.flatnonzero((Ratio_Active_Inactive_rec < Thres_Ratio) *
                                     (H[0, :] > H[1, :])
                                     )
    
Neurons_PhasicGO    = np.flatnonzero((Ratio_Active_Inactive_rec < Thres_Ratio) *
                                     (H[0, :] < H[1, :])
                                     )
    

#%% Plot components (Basis image) in 2D
plt.close("all")

Classe2Plot = [Neurons_Tonic, Neurons_PhasicNoGO, Neurons_PhasicGO]
nClasses    = len(Classe2Plot)


# Set Save dir
SaveDir = r"O:\OneDrive - MRC Laboratory of Molecular Biology\General - VGAT Paper\Figures\Paper\Figure1" + \
            r"\Version2\Maps_PhasicGONOGO_Tonic"
SaveFig = False


# Figure Params
Fig_Width   = .8 # in cm
Fig_Height  = .7 # in cm


# Axis Params
Axis_left   = -.05
Axis_bottom = -.05
Axis_width  = 1.1
Axis_height = 1.1
Axis_LineWidth  = 0

FontName = "Arial"
FontSize = 6

# Heat map params
Ticks_LineWidth = .2
Ticks_Length    = 3
YTicks          = [0, .5]


XLabel = "Time (s)"

MarkerSize = 1
    
MinFraction = 0
MaxFraction = 1


Norm            = "OneSlope"
#
Color_Tonic         = (0, 0.68, 0.94)
Color_Phasic_NoGo   =  "#fdae6b" # "#c994c7"
Color_Phasic_Go     =  "#e6550d" #"#dd1c77"


Color_Cpmnts = [Color_Tonic, Color_Phasic_NoGo, Color_Phasic_Go]
# Color_Cpmnts = ["r", "b", "g"]


ColorBar        = False


# Create Save directory
if SaveFig:
    CreateDestinationFolder(Destination_Path = SaveDir)

    
   

SaveName = rec + "_Map3d.png"


Xlim = [np.min(DB_Recs[rec]["ROIs_RowColDepth_um"][Responding_Neurons][:, 0]),
        np.max(DB_Recs[rec]["ROIs_RowColDepth_um"][Responding_Neurons][:, 0])
        ]

Ylim = [np.min(DB_Recs[rec]["ROIs_RowColDepth_um"][Responding_Neurons][:, 1]),
        np.max(DB_Recs[rec]["ROIs_RowColDepth_um"][Responding_Neurons][:, 1])
        ]

Zlim = [-np.max(DB_Recs[rec]["ROIs_RowColDepth_um"][Responding_Neurons][:, 2]),
        -np.min(DB_Recs[rec]["ROIs_RowColDepth_um"][Responding_Neurons][:, 2])
        ]

for Neuron_class in range(nClasses):
    
    
    #% Plot components (Basis image) in 3D individual
    fig = plt.figure(dpi = 300, 
                        figsize = [Fig_Width, Fig_Height],
                        label = rec
                        )
    

    ax = plt.axes([Axis_left, Axis_bottom, Axis_width, Axis_height],
                  projection='3d')
    
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", [[1, 1, 1], Color_Cpmnts[Neuron_class]])
    

    ax = PlotModel3d(np.ones(len(Classe2Plot[Neuron_class])), 
                    DB_Recs[rec]["ROIs_RowColDepth_um"][Responding_Neurons[Classe2Plot[Neuron_class]], :],
                    ax = ax,
                    MarkerSize = MarkerSize,
                    Norm = Norm, 
                    cmap = cmap, 
                    AddColorBar = ColorBar, 
                    MinFraction = MinFraction, 
                    MaxFraction = MaxFraction,
                    XTicks = [0, 500], XTickLabels = ["", ""], XLabel = '',
                    YTicks = [0, 500], YTickLabels = ["", ""], YLabel = '',
                    ZTicks = [0, -100], ZTickLabels = ["", ""], ZLabel = ''
                    )
    
    
    ax.set_xlim(Xlim)
    ax.set_ylim(Ylim)
    ax.set_zlim(Zlim)
    
    # ax.view_init(elev=30, azim=-25, roll=0)
    ax.view_init(elev=25, azim=-25, roll=0)

    # ax.axis("off")
        
    # Save figure
    if SaveFig:
        SaveName = str(Neuron_class) + ".svg"
        plt.savefig(os.path.join(SaveDir, SaveName), dpi="figure")


#%% Plot components (Basis image) in 2D superimposed


#% Plot components (Basis image) in 3D individual
fig = plt.figure(dpi = 300, 
                    figsize = [Fig_Width, Fig_Height],
                    label = rec
                    )


ax = plt.axes([Axis_left, Axis_bottom, Axis_width, Axis_height],
              projection='3d')

SaveName = rec + "_Map3d_Superimposed.svg"

for Neuron_class in range(nClasses):
    
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", [[1, 1, 1], Color_Cpmnts[Neuron_class]])
    
    ax = PlotModel3d(np.ones(len(Classe2Plot[Neuron_class])), 
                    DB_Recs[rec]["ROIs_RowColDepth_um"][Responding_Neurons[Classe2Plot[Neuron_class]], :],
                    ax = ax,
                    MarkerSize = MarkerSize,
                    Norm = Norm, 
                    cmap = cmap, 
                    AddColorBar = ColorBar, 
                    MinFraction = MinFraction, 
                    MaxFraction = MaxFraction,
                    XTicks = [0, 500], XTickLabels = ["", ""], XLabel = '',
                    YTicks = [0, 500], YTickLabels = ["", ""], YLabel = '',
                    ZTicks = [0, -100], ZTickLabels = ["", ""], ZLabel = ''
                    )

    
    ax.set_xlim(Xlim)
    ax.set_ylim(Ylim)
    ax.set_zlim(Zlim)
    
    # ax.view_init(elev=30, azim=-25, roll=0)
    ax.view_init(elev=25, azim=-25, roll=0)
    
# Save figure
if SaveFig:
    plt.savefig(os.path.join(SaveDir, SaveName), dpi="figure")

