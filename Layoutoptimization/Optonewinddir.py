#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 18:00:25 2022

@author: youssef
"""

# import os
import os

os.environ['SNOPT_LICENSE'] = '/home/youssef/snopt/snopt7.lic'
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from floris.tools import FlorisInterface
from floris.tools.optimization.layout_optimization.layout_optimization_pyoptsparse_onedir import LayoutOptimizationPyOptSparseOneDir
import pdb
from itertools import product, permutations
import matplotlib.pyplot as plt
import time
# from matplotlib.pyplot import cm
# from matplotlib.lines import Line2D
# from matplotlib.animation import FuncAnimation
def SNOPTlayoutoptimizationOneDir(fi,layout0,wd, ws, freq,optOptions,plot=1):
    fi.reinitialize(layout=layout0,
                    wind_directions=wd,
                    wind_speeds=ws
                    )
    layout_opt = LayoutOptimizationPyOptSparseOneDir(fi, boundaries, freq=freq,solver="SNOPT",
                                               optOptions=optOptions,storeHistory=None)  
    layout_opt.optimize()
    xopt, yopt=layout_opt.get_optimized_locs()
    # pdb.set_trace()
    layout=(xopt,yopt)
    if plot==1:
        layout_opt.plot_layout_opt_results()
        plt.show()
    return layout

def sort_boundaries(boundaries):
    bnd=pd.DataFrame(boundaries, columns=['x','y'])
    bnd.mean()
    A=np.degrees(np.arctan2(bnd.y-bnd.mean()[1],bnd.x-bnd.mean()[0]))
    A %= 360
    A.sort_values()
    boundaries=bnd.reindex(A.sort_values().index).values.tolist()
    return boundaries


def ConcentricCirclesLayout(N_circles,spc,D):
    spc=spc*D
    y_coordinates=np.array(0)
    x_coordinates=np.array(0)
    for i in range(N_circles):
        i=i+1
        N_turbines=math.floor(2*np.pi*i)
        angles=np.arange(0,2*np.pi,2*np.pi/N_turbines)
        x_coordinates=np.append(x_coordinates, i*spc*np.cos(angles))
        y_coordinates=np.append(y_coordinates, i*spc*np.sin(angles))
    x_coordinates=np.round(x_coordinates)
    y_coordinates=np.round(y_coordinates)
    layout= (x_coordinates.tolist(), y_coordinates.tolist())
        
    return layout

def myround(x, base=5):
    return base * round(x/base)

############## Inputs  ############## 
mindist = 6
spc = 6 # maximum spacing at the beginning of the simulation for boundary creation
shape ='square'
no_circles_squares =1 
no_turbines = 9
file_number = f'{mindist}{mindist}{no_turbines}{no_turbines}{no_turbines}'
windrose = 'iea'

Ti = 0.06

############################ 

path =fr"../Results_{windrose}/Results_{mindist}_{spc}"
angle=1
opt_baseline_file ='/SNOPTlayoutnew_'+shape+'_'+str(no_turbines)+'T_iea.csv'
sum_print='SNOPTlayoutnew_'+shape+'_'+str(no_turbines)+'T_iea_TI_'+str(Ti)+'.out'
opt_dyn_layout_file='/OptimumDynamicfarmlayout_'+shape+'_'+str(no_turbines)+'T_iea_TI_'+str(Ti)+'.csv'
opt_dyn_layout_file_xaxis ='/OptimumDynamicfarmlayout_xaxis_'+shape+'_'+str(no_turbines)+'T_iea_TI_'+str(Ti)+'.csv'
opt_dyn_layout_file_yaxis ='/OptimumDynamicfarmlayout_yaxis_'+shape+'_'+str(no_turbines)+'T_iea_TI_'+str(Ti)+'.csv'
opt_aep_sum_file ='/AEP_summary_per_winddir_'+shape+'_'+str(no_turbines)+'T_iea_TI_'+str(Ti)+'.csv'



# Initialize the FLORIS interface fi
fi = FlorisInterface('../Inputs/gch15MW.yaml')
# pdb.set_trace()
D = fi.floris.farm.rotor_diameters_sorted[0][0][0]
fi.floris.flow_field.turbulence_intensity=Ti

############## Circle  ############## 
if shape=='circle':
    N_circles=no_circles_squares
    spacing=spc*D
    spc=spc
    angles=np.arange(0,360+angle,angle)
    boundaries_x=np.round(((N_circles*spacing)+0.01*D)*np.cos(np.radians(angles)))
    boundaries_y=np.round(((N_circles*spacing)+0.01*D)*np.sin(np.radians(angles)))
    boundaries = [[x,y] for x, y in zip(boundaries_x, boundaries_y)]
    # boundaries=sort_boundaries(boundaries)
    layout0=ConcentricCirclesLayout(N_circles,spc,D)    


############## Square  ############## 
## D=240
if shape=='square':
    N_squares=no_circles_squares
    spacing=spc*D
    y_coordinates=np.tile(np.arange(-spacing*N_squares,spacing*N_squares+spacing,spacing),1+N_squares*2)
    x_coordinates=np.repeat(np.arange(-spacing*N_squares,spacing*N_squares+spacing,spacing),1+N_squares*2)
    x_coordinates=np.round(x_coordinates)
    y_coordinates=np.round(y_coordinates)
    layout0= (x_coordinates.tolist(), y_coordinates.tolist())
    
    boundaries_x=np.round((min(layout0[0])-0.01*D,min(layout0[0])-0.01*D,
                            max(layout0[0])+0.01*D,max(layout0[0])+0.01*D,
                            min(layout0[0])-0.01*D),decimals=2)
    boundaries_y=np.round((min(layout0[0])-0.01*D,max(layout0[0])+0.01*D,
                            max(layout0[0])+0.01*D,min(layout0[0])-0.01*D,
                            min(layout0[0])-0.01*D),decimals=2)
    boundaries = [[x,y] for x, y in zip(boundaries_x, boundaries_y)]
    boundaries = pd.DataFrame(boundaries, columns=['x','y'])
    boundaries=boundaries.values.tolist()

################Horns Rev############  

if shape=='HornsRevall3':
    
    HornsRev_layout= pd.read_csv('../Inputs/HornsRev_layout.csv')
    x = HornsRev_layout['x']
    y = HornsRev_layout['y']
    layout0= (x.tolist(), y.tolist())
    
    # boundaries_x=np.round((-6843-0.01*D,min(layout0[0])-0.01*D,
    #                         6843,max(layout0[0])+0.01*D,
    #                         -6843-0.01*D),decimals=2)
    # boundaries_y=np.round((min(layout0[1])-0.01*D,max(layout0[1])+0.01*D,
    #                         max(layout0[1])+0.01*D,min(layout0[1])-0.01*D,
    #                         min(layout0[1])-0.01*D),decimals=2)
    boundaries_x=np.round((-6843,min(layout0[0]),
                            6843,max(layout0[0]),
                            -6843),decimals=2)
    boundaries_y=np.round((min(layout0[1]),max(layout0[1]),
                            max(layout0[1]),min(layout0[1]),
                            min(layout0[1])),decimals=2)
    boundaries = [[x,y] for x, y in zip(boundaries_x, boundaries_y)]
    boundaries = pd.DataFrame(boundaries, columns=['x','y'])
    boundaries=boundaries.values.tolist()
    
# fig = plt.figure(figsize=(5,5),dpi=96)
# ax = fig.add_subplot()
# ax.plot(boundaries_x,boundaries_y)
# ax.scatter(x,y)
# plt.show()
############################  

if shape=='HornsRevall3':
    layout=layout0
    del layout0
else:  
    layoutdf=pd.read_csv(path+ '/BaselineOptimization' +opt_baseline_file)
    layout=(np.array(layoutdf.x).tolist(),
            np.array(layoutdf.y).tolist())

# Windrose parameters
if windrose =='iea':
    print('############################ iea windrose ############################')
    wd=[0., 22.5, 45., 67.5, 90., 112.5, 135., 157.5, 180., 202.5, 225., 247.5, 270., 292.5, 315., 337.5]
    ws=[10, 10,  10,  10,    10,   10,    10,   10,    10,    10,    10,   10,   10,    10,   10,    10]
    freq_d= [.025,  .024,  .029,  .036,.063,  .065,  .100,  .122,.063,  .038,  .039,  .083, .213,  .046,  .032,  .022]
    ########################################################
if windrose =='alpha':
    print(' ########################### alpha ventus windrose ############################')
    wd=[0., 22.5, 45., 67.5, 90., 112.5, 135., 157.5, 180., 202.5, 225., 247.5, 270., 292.5, 315., 337.5]
    ws=[10, 10,  10,  10,    10,   10,    10,   10,    10,    10,    10,   10,   10,    10,   10,    10]
    freq_d=[0.0313,0.0402,0.0375, 0.0568,0.0558,0.0608, 0.0424,0.0564,0.0555, 0.1114,0.0932,0.1114, 0.0722,
          0.0743,0.0500, 0.0508]
    
if windrose =='HornsRev4':
    print(' ########################### HornsRev windrose ############################')
    f2 = [3.597152, 3.948682, 5.167395, 7.000154, 8.364547, 6.43485,
                  8.643194, 11.77051, 15.15757, 14.73792, 10.01205, 5.165975]
    # freq_d = np.array(f) / np.sum(f)
    # wd=[0, 30 ,60, 90, 120, 150, 180.,  210, 240, 270, 300 , 330]
    # ws=[10, 10,  10,  10,    10,   10,    10,   10,    10,    10,    10,   10 ]
    
    f = [7.000154, 14.73792]
    freq_d = np.array(f) / np.sum(f2)
    wd=[90, 270]
    ws=[10, 10]

if windrose =='HornsRevall3':
    print(' ########################### HornsRev windrose all############################')
    f2 = [3.597152, 3.948682, 5.167395, 7.000154, 8.364547, 6.43485,
                  8.643194, 11.77051, 15.15757, 14.73792, 10.01205, 5.165975]
    # freq_d = np.array(f) / np.sum(f)
    # wd=[0, 30 ,60, 90, 120, 150, 180.,  210, 240, 270, 300 , 330]
    # ws=[10, 10,  10,  10,    10,   10,    10,   10,    10,    10,    10,   10 ]
    
    f = [3.597152, 3.948682, 5.167395,  8.364547, 6.43485,
                  8.643194, 11.77051, 15.15757,  10.01205, 5.165975]
    freq_d = np.array(f) / np.sum(f2)
    wd=[0, 30 ,60,  120, 150, 180.,  210, 240,  300 , 330]
    ws=[10, 10,  10,  10,    10,   10,    10,   10,    10,    10]


    # #######################################################

Windrose=pd.DataFrame()
Windrose['wd']=wd
Windrose['ws']=ws
Windrose['freq']=freq_d
Windrose['MoorPy_angles']=90-Windrose.wd+180
Windrose['MoorPy_angles']%= 360


layoutx=pd.DataFrame(np.tile(np.array(layout[0]),(len(Windrose['MoorPy_angles']),1)).transpose(),columns=Windrose['MoorPy_angles'])
layouty=pd.DataFrame(np.tile(np.array(layout[1]),(len(Windrose['MoorPy_angles']),1)).transpose(),columns=Windrose['MoorPy_angles'])
df=pd.DataFrame(np.tile(-Windrose['MoorPy_angles'],(len(layoutx),1)), columns=layoutx.columns)
layoutx_new=round(layoutx*np.cos(np.deg2rad(-df))+layouty*np.sin(np.deg2rad(-df)),2)
layouty_new=round(layouty*np.cos(np.deg2rad(-df))-layoutx*np.sin(np.deg2rad(-df)),2)

boundaries_x=pd.DataFrame(np.tile(np.array(boundaries_x),(len(Windrose['MoorPy_angles']),1)).transpose(),columns=Windrose['MoorPy_angles'])
boundaries_y=pd.DataFrame(np.tile(np.array(boundaries_y),(len(Windrose['MoorPy_angles']),1)).transpose(),columns=Windrose['MoorPy_angles'])
df=pd.DataFrame(np.tile(-Windrose['MoorPy_angles'],(len(boundaries_x),1)), columns=layoutx.columns)
boundaries_x_new=round(boundaries_x*np.cos(np.deg2rad(-df))+boundaries_y*np.sin(np.deg2rad(-df)),2)
boundaries_y_new=round(boundaries_y*np.cos(np.deg2rad(-df))-boundaries_x*np.sin(np.deg2rad(-df)),2)


AEP_summary=pd.DataFrame(Windrose['MoorPy_angles'])
AEP_summary['MoorPy_angles_rnd']=myround(Windrose['MoorPy_angles'],base=5)
Floris_perp=pd.DataFrame()
Floris_X=pd.DataFrame()
Floris_Y=pd.DataFrame()
for i in range(len(AEP_summary)):
    # i=7
    layout0=(np.array(layoutx_new.iloc[:,i]).tolist(),
            np.array( layouty_new.iloc[:,i]).tolist())   
    boundaries = [[x,y] for x, y in zip(np.array(boundaries_x_new.iloc[:,i]).tolist(), 
                                        np.array( boundaries_y_new.iloc[:,i]).tolist())]
    # pdb.set_trace()
    wd=[270]
    ws=[10] 
    freq= [Windrose['freq'][i]]
    freq=np.transpose(np.array([freq]))
    # fi.reinitialize_flow_field(layout_array=layout0)  
    fi.reinitialize(layout=layout0,wind_directions=wd,wind_speeds=ws)
    AEP_initial = fi.get_farm_AEP(freq=freq,freq_warn=0)* 1e-9
    print("=====================================================")
    print('AEP_initial='+str(AEP_initial))
    print("=====================================================")
    AEP_summary.loc[i,'AEP_initial']=AEP_initial
    optOptions={"Major feasibility tolerance": 1e-6, 
                   # "Verify level": 3, 
                  "Scale option":0 ,
                  # "Function precision": 1e-6,
                "Major optimality tolerance": 1e-5,
                # "Major optimality tolerance": 1e-7,  #what i will use for PhD isA
                # "Major optimality tolerance": 1e-6,  #what i will use for PhD isA
                # "Derivative level":1, 
                "iPrint": int(file_number)-1,
                "iSumm":  int(file_number),
                'Print file':'print_'+sum_print,
                'Summary file':'sum_'+sum_print,
                # "Linesearch tolerance": 0.01
                }
    p=0
    # if i ==1:
    #     p=1
    layout=SNOPTlayoutoptimizationOneDir(fi,layout0,wd, ws, freq,optOptions,plot=p)
    fi.reinitialize(layout=layout)
    AEP_new= fi.get_farm_AEP(freq=freq,freq_warn=0)* 1e-9
    print("=====================================================")
    print('AEP_new='+str(AEP_new))
    print("=====================================================")
    AEP_summary.loc[i,'AEP_new']=AEP_new
    Floris_perp.loc[:,AEP_summary['MoorPy_angles_rnd'][i]]=np.array(layout[1])-np.array(layout0[1])
    # pdb.set_trace()
    ang=AEP_summary['MoorPy_angles'][i]
    Floris_X.loc[:,AEP_summary['MoorPy_angles_rnd'][i]]=np.array(layout[0])*np.cos(np.deg2rad(ang))-np.array(layout[1])*np.sin(np.deg2rad(ang))
    Floris_Y.loc[:,AEP_summary['MoorPy_angles_rnd'][i]]=np.array(layout[1])*np.cos(np.deg2rad(ang))+np.array(layout[0])*np.sin(np.deg2rad(ang))

###############
AEP_summary['AEP_gain']=(AEP_summary.AEP_new-AEP_summary.AEP_initial)/np.sum(AEP_summary.AEP_initial)*100
AEP_summary.sort_values(by=['AEP_gain'],ignore_index=True,inplace=True,ascending=False)
Floris_perp=Floris_perp.sort_index(axis=1)
Floris_X= Floris_X.sort_index(axis=1)
Floris_Y= Floris_Y.sort_index(axis=1)

AEP_summary.to_csv(path+ '/SingleWIndDir_Opt'+ opt_aep_sum_file , index=False)
Floris_perp.to_csv(path+ '/SingleWIndDir_Opt'+opt_dyn_layout_file, index=False)
Floris_X.to_csv(path+ '/SingleWIndDir_Opt'+opt_dyn_layout_file_xaxis, index=False)
Floris_Y.to_csv(path+ '/SingleWIndDir_Opt'+opt_dyn_layout_file_yaxis, index=False)
########### Test to verify the coordinate transformation

# plt.axis("equal")
# plt.grid()
# plt.tick_params(which="both", labelsize=fontsize)
# plt.show()