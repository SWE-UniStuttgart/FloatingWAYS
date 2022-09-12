#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 27 13:00:07 2022

@author: youssef
"""

import os
os.environ['SNOPT_LICENSE'] = '/home/youssef/snopt/snopt7.lic'
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import floris.tools as wfct
import floris.tools.optimization.pyoptsparse as opt
import pdb
from itertools import product, permutations
import matplotlib.pyplot as plt


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


def SNOPTlayoutoptimization(fi,layout0,wd, ws, freq,plot):
    model = opt.layout.Layout(fi, boundaries, wdir=wd, wspd=ws, wfreq=freq)
    # optOptions={"Major feasibility tolerance": 1e-6, "Verify level": 3, "Scale option":2 ,"Major optimality tolerance": 5e-5}
    tmp = opt.optimization.Optimization(model=model, solver="SNOPT")    
    sol = tmp.optimize()
    if plot==1: 
        model.plot_layout_opt_results(sol)
        plt.show()    
        
    layout=(sol.getDVs()["x"].tolist(),sol.getDVs()["y"].tolist())
    
    return layout

def savelayout(layout,path,filename):
    layoutdf=pd.DataFrame(layout[0],columns=['x'])
    layoutdf['y']=layout[1]
    layoutdf.to_csv(path+filename, index=False)




# Initialize the FLORIS interface fi
fi = wfct.floris_interface.FlorisInterface("FLORIS_15MW.json")
D = fi.floris.farm.turbines[0].rotor_diameter

N_circles=2
spacing=5*D
spc=5
angles=np.arange(0,360,1)
boundaries_x=np.round(N_circles*spacing*np.cos(np.radians(angles)))
boundaries_y=np.round(N_circles*spacing*np.sin(np.radians(angles)))
boundaries = [[x,y] for x, y in zip(boundaries_x, boundaries_y)]
boundaries=sort_boundaries(boundaries)
layout0=ConcentricCirclesLayout(N_circles,spc,D)    

wd=[0., 22.5, 45., 67.5, 90., 112.5, 135., 157.5, 180., 202.5, 225., 247.5, 270., 292.5, 315., 337.5]
ws=[10, 10,  10,  10,    10,   10,    10,   10,    10,    10,    10,   10,   10,    10,   10,    10]
freq= [.025,  .024,  .029,  .036,.063,  .065,  .100,  .122,.063,  .038,  .039,  .083, .213,  .046,  .032,  .022]

fi.reinitialize_flow_field(layout_array=layout0)
AEP_initial=fi.get_farm_AEP(np.array(wd), np.array(ws), np.array(freq)) * 1e-9
print("=====================================================")
print('AEP_initial='+str(AEP_initial))
print("=====================================================")

layout=SNOPTlayoutoptimization(fi,layout0,wd, ws, freq,1)

fi.reinitialize_flow_field(layout_array=layout)
AEP_optimized=fi.get_farm_AEP(np.array(wd), np.array(ws), np.array(freq)) * 1e-9
print("=====================================================")
print('AEP_current='+str(AEP_optimized))
print("=====================================================")

savelayout(layout,'','SNOPTlayoutnew.csv')