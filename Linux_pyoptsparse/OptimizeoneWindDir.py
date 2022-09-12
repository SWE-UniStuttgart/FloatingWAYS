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
import floris.tools as wfct
import floris.tools.optimization.pyoptsparse as opt
import pdb
from itertools import product, permutations
import matplotlib.pyplot as plt
# from matplotlib.pyplot import cm
# from matplotlib.lines import Line2D
# from matplotlib.animation import FuncAnimation

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



fi = wfct.floris_interface.FlorisInterface(
    "FLORIS_15MW.json"
)
D = fi.floris.farm.turbines[0].rotor_diameter

N_circles=2
spacing=5*D
spc=5
angles=np.arange(0,360,1)
boundaries_x=np.round(N_circles*spacing*np.cos(np.radians(angles)))
boundaries_y=np.round(N_circles*spacing*np.sin(np.radians(angles)))
boundaries = [[x,y] for x, y in zip(boundaries_x, boundaries_y)]
boundaries=sort_boundaries(boundaries)
# layout=ConcentricCirclesLayout(N_circles,spc,D)


layout=(np.array([ 1323.01,  1574.75,   358.74,  -249.84 ,
        -820.66, -2275.84, -1700.47,  2279.15,
        2013.42,   232.56,  -522.56, -1290.64,
       -2004.19, -2339.12, -1064.76, -1505.33,
         -17.96,   883.87,  1953.94]).tolist(),
        np.array([ -460.14,   647.74,   929.29, -2387.11,
         160.84,  -761.93,  1406.29,   186.07,
        1306.33,  1733.7,  2252.14,  2023.28,
         423.5 ,  -161.59, -1052.92, -1801.91,
        -851.32, -2146.11, -1394.09 ]).tolist())

wd=[0., 22.5, 45., 67.5, 90., 112.5, 135., 157.5, 180., 202.5, 225., 247.5, 270., 292.5, 315., 337.5]
ws=[10, 10,  10,  10,    10,   10,    10,   10,    10,    10,    10,   10,   10,    10,   10,    10]
freq= [.025,  .024,  .029,  .036,.063,  .065,  .100,  .122, .063,  .038,  .039,  .083, .213,  .046,  .032,  .022]



Windrose=pd.DataFrame()
Windrose['wd']=wd
Windrose['ws']=ws
Windrose['freq']=freq
Windrose['MoorPy_angles']=90-Windrose.wd+180
Windrose['MoorPy_angles']%= 360


layoutx=pd.DataFrame(np.tile(np.array(layout[0]),(len(Windrose['MoorPy_angles']),1)).transpose(),columns=Windrose['MoorPy_angles'])
layouty=pd.DataFrame(np.tile(np.array(layout[1]),(len(Windrose['MoorPy_angles']),1)).transpose(),columns=Windrose['MoorPy_angles'])
df=pd.DataFrame(np.tile(-Windrose['MoorPy_angles'],(len(layoutx),1)), columns=layoutx.columns)
layoutx_new=round(layoutx*np.cos(np.deg2rad(-df))+layouty*np.sin(np.deg2rad(-df)),2)
layouty_new=round(layouty*np.cos(np.deg2rad(-df))-layoutx*np.sin(np.deg2rad(-df)),2)

# try_x=layoutx_new*np.cos(np.deg2rad(-df))-layouty_new*np.sin(np.deg2rad(-df))
# try_y=layouty_new*np.cos(np.deg2rad(-df))+layoutx_new*np.sin(np.deg2rad(-df))

AEP_summary=pd.DataFrame(Windrose['MoorPy_angles'])
AEP_summary['MoorPy_angles_rnd']=myround(Windrose['MoorPy_angles'],base=5)
Floris_perp=pd.DataFrame()
Floris_X=pd.DataFrame()
Floris_Y=pd.DataFrame()
for i in range(len(AEP_summary)):
    # i=7
    layout0=(np.array(layoutx_new.iloc[:,i]).tolist(),
            np.array( layouty_new.iloc[:,i]).tolist())   
    # pdb.set_trace()
    wd=[270] #This is the zero wind direction according to floris
    ws=[10] 
    freq= [Windrose['freq'][i]]
    fi.reinitialize_flow_field(layout_array=layout0)    
    AEP_initial=fi.get_farm_AEP(np.array(wd), np.array(ws), np.array(freq)) * 1e-9
    print("=====================================================")
    print('AEP_initial='+str(AEP_initial))
    print("=====================================================")
    AEP_summary.loc[i,'AEP_initial']=AEP_initial
    model = opt.layoutOneWindDir.Layout(fi, boundaries, wdir=wd, wspd=ws, wfreq=freq)
    optOptions={"Major feasibility tolerance": 1e-6, "Verify level": 3, "Scale option":2 ,"Major optimality tolerance": 5e-5}
    tmp = opt.optimization.Optimization(model=model, solver="SNOPT")    
    sol = tmp.optimize()
    model.plot_layout_opt_results(sol)
    plt.show()    
    layout=(np.array(layoutx_new.iloc[:,i]).tolist(),sol.getDVs()["y"].tolist())
    fi.reinitialize_flow_field(layout_array=layout)
    AEP_new=fi.get_farm_AEP(np.array(wd), np.array(ws), np.array(freq)) * 1e-9
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

AEP_summary.to_csv('AEP_summary_per_winddir.csv', index=False)
Floris_perp.to_csv('OptimumDynamicfarmlayout.csv', index=False)
Floris_X.to_csv('OptimumDynamicfarmlayout_xaxis.csv', index=False)
Floris_Y.to_csv('OptimumDynamicfarmlayout_yaxis.csv', index=False)
########### Test to verify the coordinate transformation
# fi.reinitialize_flow_field(layout_array=layout)
# AEP=pd.DataFrame()

# for i in range(16):
#     layout=(np.array(layoutx_new.iloc[:,i]).tolist(),
#             np.array( layouty_new.iloc[:,i]).tolist())
#     fi.reinitialize_flow_field(layout_array=layout)
#     wd=[270]
#     ws=[10]
#     freq= [Windrose['freq'][i]]
#     AEP_initial=fi.get_farm_AEP(np.array(wd), np.array(ws), np.array(freq)) * 1e-9
#     AEP.loc[i,1]=AEP_initial



# plt.figure(figsize=(9, 6))
# fontsize = 16
# plt.plot(layout0[0], layout0[1], "ob")
# plt.plot(layout[0], layout[1], "or")
# # plt.title('Layout Optimization Results', fontsize=fontsize)
# plt.xlabel("x (m)", fontsize=fontsize)
# plt.ylabel("y (m)", fontsize=fontsize)
# plt.axis("equal")
# plt.grid()
# plt.tick_params(which="both", labelsize=fontsize)
# plt.show()