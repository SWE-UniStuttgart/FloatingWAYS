#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 11 07:52:00 2023

@author: youssef
"""


import os
os.environ['SNOPT_LICENSE'] = '/home/youssef/snopt/snopt7.lic'
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from floris.tools import FlorisInterface
from floris.tools.optimization.layout_optimization.layout_optimization_pyoptsparse import LayoutOptimizationPyOptSparse
import pdb
# from floris.tools import visualize_cut_plane #, plot_turbines_with_fi
# import floris.tools as wfct
# import shapely
from shapely.geometry import Polygon, LineString, Point


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


def SNOPTlayoutoptimization(fi,layout0,wd, ws, freq,optOptions,minimum_D):
    fi.reinitialize(layout=layout0,
                    wind_directions=wd,
                    wind_speeds=ws)
    layout_opt = LayoutOptimizationPyOptSparse(fi, boundaries, freq=freq,solver="SNOPT",
                                               optOptions=optOptions,storeHistory=None,min_dist=minimum_D)  
    layout_opt.optimize()
    xopt, yopt=layout_opt.get_optimized_locs()
    layout=(xopt,yopt)
    
    layout_opt.plot_layout_opt_results()
    plt.show()
    return layout

def savelayout(layout,path,filename):
    layoutdf=pd.DataFrame(layout[0],columns=['x'])
    layoutdf['y']=layout[1]
    layoutdf.to_csv(path+filename, index=False)


############## Inputs  ############## 
mindist = 3
spc = 5 # maximum spacing at the beginning of the simulation for boundary creation
shape ='square'
no_circles_squares =1 
no_turbines = 9 
file_number = f'{mindist}{mindist}{mindist}{no_turbines}{no_turbines}{no_turbines}'
windrose = 'iea'

Ti = 0.06

############################ 

path =fr"../Results_{windrose}/Results_{mindist}_{spc}"
angle=1


opt_baseline_file =f'/SNOPTlayoutnew_{shape}_{no_turbines}T_{windrose}.csv'
sum_print=f'SNOPTlayoutnew_{shape}_{no_turbines}T_{windrose}.out'


# Initialize the FLORIS interface fi
fi = FlorisInterface('../Inputs/gch15MW.yaml')
D = fi.floris.farm.rotor_diameters_sorted[0][0][0]
fi.floris.flow_field.turbulence_intensity=Ti

############## Circle  ############## 
if shape=='circle':
    N_circles=no_circles_squares
    spacing=spc*D
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
    spc=spc*D
    y_coordinates=np.tile(np.arange(-spc*N_squares,spc*N_squares+spc,spc),1+N_squares*2)
    x_coordinates=np.repeat(np.arange(-spc*N_squares,spc*N_squares+spc,spc),1+N_squares*2)
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
    boundaries = boundaries.values.tolist()


boundarypolygon = Polygon(boundaries)
boundaryline = LineString(boundaries)
boundarycon = np.zeros(len(layout0[0]))
for i in range(len(layout0[0])):
    loc = Point(layout0[0][i], layout0[1][i])
    boundarycon[i] = loc.distance(boundaryline)
    if boundarypolygon.contains(loc)==True:
        boundarycon[i] *= -1.0
        
if not all(boundarycon<0):
    print('not all turbines lie within the boundaries at initialization')
###########################  

############################ iea task ############################
if windrose =='iea':
    print('############################ iea windrose ############################')
    wd=[0., 22.5, 45., 67.5, 90., 112.5, 135., 157.5, 180., 202.5, 225., 247.5, 270., 292.5, 315., 337.5]
    ws=[10]
    freq_d= [.025,  .024,  .029,  .036,.063,  .065,  .100,  .122,.063,  .038,  .039,  .083, .213,  .046,  .032,  .022]
    ########################################################
if windrose =='alpha':
    print(' ########################### alpha ventus windrose ############################')
    wd=[0., 22.5, 45., 67.5, 90., 112.5, 135., 157.5, 180., 202.5, 225., 247.5, 270., 292.5, 315., 337.5]
    ws=[10]
    freq_d=[0.0313,0.0402,0.0375, 0.0568,0.0558,0.0608, 0.0424,0.0564,0.0555, 0.1114,0.0932,0.1114, 0.0722,
          0.0743,0.0500, 0.0508]
    # #######################################################
    
freq=np.transpose(np.array([freq_d]))
fi.reinitialize(layout=layout0,
                wind_directions=wd,
                wind_speeds=ws)

AEP_initial = fi.get_farm_AEP(freq=freq)* 1e-9
print("=====================================================")
print('AEP_initial='+str(AEP_initial))
print("=====================================================")

optOptions={"Major feasibility tolerance": 1e-6, 
               # "Verify level": 3, 
              "Scale option":0 ,
              # "Function precision": 1e-6,
            # "Major optimality tolerance": 5e-5,
            "Major optimality tolerance": 1e-5,
            # "Derivative level":1, 
            "iPrint": int(file_number)-1,
            "iSumm":  int(file_number),
            'Print file':'print_'+sum_print,
            'Summary file':'sum_'+sum_print,
            # "Linesearch tolerance": 0.01
            }
layout=SNOPTlayoutoptimization(fi,layout0,wd, ws, freq,optOptions,minimum_D=mindist*D)
fi.reinitialize(layout=layout)
AEP_optimized= fi.get_farm_AEP(freq=freq)* 1e-9
print("=====================================================")
print('AEP_current='+str(AEP_optimized))
print("=====================================================")


savelayout(layout, path+ '/BaselineOptimization' , opt_baseline_file)