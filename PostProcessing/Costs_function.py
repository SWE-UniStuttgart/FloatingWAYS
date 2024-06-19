#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 14:26:37 2024

@author: youssef
"""



import moorpy as mp
import os
import numpy as np
import pandas as pd
import pdb
from celluloid import Camera
import matplotlib.pyplot as plt
import itertools
from shapely.geometry import LineString
import string

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size" :10,
    "axes.titlesize":10,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
})

def LCOE(d,L):
    MBL = 30.2*d**2*(44-80*d)*10**6
    Cost = ((0.0591*MBL*10**-3-87.6)*L) + (10.198 *MBL*10**-3)
    return sum(Cost)
########## Inputs Start ##########

########### Simulation input
mindist = 4 #this is the small number
spc = 6 # maximum spacing at the beginning of the simulation for boundary creation
shape ='square'
no_circles_squares =1 
no_turbines = 9 
windrose = 'iea'

tolerance= np.array([ 5,10,20,
            30,40,50,
            # 5,10,20,
            60,70,80,90,100,110,120,
            130,
            140,150,160,170,180,190,200,
            # 210,220,230,240
            ])
filterMS = 'designheading' #designonly or designheading  
MooringDatabase = 'L04' # L04 or L03. L03 means line length 03,05,07, anchor radius 2.5, 3D. L04 line length 04,06,08 anc radius 2.5, 3.5D
WindSpeed=10

######### Baseline cost ###################
d=np.array([0.1902,0.1902,0.1902])
L=np.array([623.3,623.3,623.3])
Cost_base=9*LCOE(d,L)
############################ 
path0=fr"../Results_{windrose}"
path1 =fr"../Results_{windrose}/Results_{mindist}_{spc}"
path =path1 + r"/FinalMooringWFDesign"
MooringPermMatrix= pd.read_csv(path0+f'/MooringDatabase_{MooringDatabase}/WindSpeed_10'+ '/Final_PermMatrix.csv') # the design parameters vof each mooring system
resultcost = pd.DataFrame()
counter=0
for j in tolerance:
    tol = j
    outputfile=f'/TBMDfinal_tol_{tol}_filter_{filterMS}_MooringDatabase_{MooringDatabase}'+'.csv'
    layout_final=pd.read_csv(path + outputfile,index_col=False)
    
    Farm_Matrix = MooringPermMatrix.loc[layout_final.Design,:].reset_index(drop=True)
    # Farm_Matrix=pd.read_csv(path+f'/Turbines_tol_{tol}_filter_{filterMS}_MooringDatabase_{MooringDatabase}/WindSpeed_{WindSpeed}/X_PermMatrix.csv')
    
    
    Cost_Matrix = pd.DataFrame()
    
    for i in Farm_Matrix.index:
        d=np.array([Farm_Matrix.D1[i], Farm_Matrix.D2[i], Farm_Matrix.D3[i]])
        L=np.array([Farm_Matrix.L1[i], Farm_Matrix.L2[i], Farm_Matrix.L3[i]])
        Cost_Matrix.loc[i,0]=LCOE(d,L)
    
    print((-Cost_Matrix.sum()[0]+Cost_base)/Cost_base*100)
    resultcost.loc[counter,'tolerance'] = tol
    resultcost.loc[counter,'cost'] = (-Cost_Matrix.sum()[0]+Cost_base)/Cost_base*100
    counter = counter + 1

import seaborn as sns
mydpi=200
width=450/72
colors = sns.color_palette()
x = np.arange(len(resultcost)) 

fig, ax = plt.subplots(figsize=(width,width*(5**.5 - 1) / 2*0.5),dpi=mydpi)
# offset = width * multiplier
ax.bar(x ,resultcost['cost'], 0.8, label=[], color = 'r')
print(ax.get_label())  
# multiplier += 1
    
handles, labels=ax.get_legend_handles_labels()
legend=dict(zip(labels,handles))
ax.grid(True)
ax.set_ylabel('Gain of MS cost [\%]')
ax.set_xlabel('tolerance [m]')
ax.set_xticks(x , ['5', '10', '20', '30', 
                    '40', '50', '60', '70',
                    '80', '90', '100', '110', '120',
                    '130','140','150','160', '170',
                    '180', '190', '200'],
              rotation=90)
keys=legend.keys()
ax.grid(True)
plt.grid(True)
plt.savefig(r"/media/sf_Shared_folder"+'/PhD_Figs'+ '/numdesign_cost.png',dpi=200,bbox_inches='tight')