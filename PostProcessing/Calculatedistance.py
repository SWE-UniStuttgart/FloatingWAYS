#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 18 13:28:00 2023

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



def rotationMatrix(x3,x2,x1):
    '''Calculates a rotation matrix based on order-z,y,x instrinsic (tait-bryan?) angles, meaning
    they are about the ROTATED axes. (rotation about z-axis would be (0,0,theta) )
    
    Parameters
    ----------
    x3, x2, x1: floats
        The angles that the rotated axes are from the nonrotated axes. Normally roll,pitch,yaw respectively. [rad]
    Returns
    -------
    R : matrix
        The rotation matrix
    '''
    # initialize the sines and cosines
    s1 = np.sin(x1) 
    c1 = np.cos(x1)
    s2 = np.sin(x2) 
    c2 = np.cos(x2)
    s3 = np.sin(x3) 
    c3 = np.cos(x3)
    
    # create the rotation matrix
    R = np.array([[ c1*c2,  c1*s2*s3-c3*s1,  s1*s3+c1*c3*s2],
                  [ c2*s1,  c1*c3+s1*s2*s3,  c3*s1*s2-c1*s3],
                  [   -s2,           c2*s3,           c2*c3]])
    
    return R

def FAero_windD(Faero,Maero,windD,rollangle=0,pitchangle=0,yawangle=0):
    # windD=np.arange(0,360,5)
    # rollangle=2
    # pitchangle=4
    # yawangle=0
    # Faero=np.array([2.149e+6,0,0,1.886e+7,2.089e+8,1.046e+6])
    # Mxdiff0=3.102e+7
    # Mydiff0=2.96e+7
    # Mzdiff0=1142e+3
    
    
    rollanglenew=rollangle*np.cos(np.deg2rad(windD))+pitchangle*np.sin(np.deg2rad(windD)) # wind direction
    pitchanglenew=pitchangle*np.cos(np.deg2rad(windD))-rollangle*np.sin(np.deg2rad(windD))

    
    
    uptilt=6
    rHub=np.array([-11.075*np.cos(np.deg2rad(uptilt))+0.2736 + 0.0439*pitchanglenew, 
                  -0.0588 - 0.048*rollanglenew, 
                  135.0005])
    

    # Faero=np.array([2.152e+6,-5.612e3,-1.499e+4])
    # Maero=np.array([1.718e7,3.872e6,2.706e6])
    
    Faero_new=Faero.copy()
    Maero_new=Maero.copy()
    ############## Works with old system ##############
    Maero_new[0]=Maero[0] - 8.5e4*pitchanglenew
    # Maero[1]=3.872e6 - 3.6e5*pitchanglenew - 3e4*yawangle + 4.3e4*rollanglenew
    Maero_new[2]=Maero[2] + 1e4*pitchanglenew  - 2e5*yawangle - 1e5*rollanglenew #-3e5*yawangle+ 1e4*pitchanglenew 
    # Faero[0] = 2.152e+6 - 1.25e4*pitchanglenew -3e3 *abs(yawangle)
    Faero_new[1] = Faero[1] - 2e3*yawangle + 9e2*pitchanglenew - 5e2*rollanglenew
    Faero_new[2] = Faero[2] + 9e2*yawangle + 2.2e3*pitchanglenew + 2e2*rollanglenew
    ############## Works with old system ##############

    Faero2=Faero_new.copy()
    Maero2=Maero_new.copy()
    Faero2[0]=Faero_new[0]*np.cos(np.deg2rad(uptilt))+Faero_new[2]*np.sin(np.deg2rad(uptilt))
    Faero2[2]=Faero_new[2]*np.cos(np.deg2rad(uptilt))-Faero_new[0]*np.sin(np.deg2rad(uptilt))


    Maero2[0]=Maero_new[0]*np.cos(np.deg2rad(uptilt))+Maero_new[2]*np.sin(np.deg2rad(uptilt))
    Maero2[2]=Maero_new[2]*np.cos(np.deg2rad(uptilt))-Maero_new[0]*np.sin(np.deg2rad(uptilt))

    
    Ang=np.array([np.deg2rad(0.0818+0.0426*rollanglenew), np.deg2rad(0.1543+0.039*pitchanglenew), np.deg2rad(0)])
    RotMat = rotationMatrix(Ang[0], -Ang[1], Ang[2]) 
    Faero3 = np.matmul(RotMat,Faero2)
    Maero3 = np.matmul(RotMat,Maero2)
    


    Ang=np.array([np.deg2rad(rollanglenew), np.deg2rad(pitchanglenew), np.deg2rad(yawangle)])
    RotMat = rotationMatrix(Ang[0], Ang[1], Ang[2]) 
    F_rot = np.matmul(RotMat,Faero3)
    M_rot = np.matmul(RotMat,Maero3)
    

    # F_rot=Faero
    # M_rot=Maero
    rRel = np.matmul(RotMat,rHub)
    moment=np.cross(rRel , F_rot)
    Mfinal= M_rot + moment
    # F_rot[0]=F_rot[0]+51000.0
    

    
    Mxdiff_case=Mfinal[0]
    Mydiff_case=Mfinal[1]
    Mzdiff_case=Mfinal[2]
    ######################## 
    Mxdiff_casenew=Mxdiff_case*np.cos(np.deg2rad(-windD))+Mydiff_case*np.sin(np.deg2rad(-windD)) # negative wind direction
    Mydiff_casenew=Mydiff_case*np.cos(np.deg2rad(-windD))-Mxdiff_case*np.sin(np.deg2rad(-windD))
    
    F=np.array([F_rot[0]*np.cos(np.deg2rad(-windD))+F_rot[1]*np.sin(np.deg2rad(-windD)),
                  F_rot[1]*np.cos(np.deg2rad(-windD))-F_rot[0]*np.sin(np.deg2rad(-windD)),
                  F_rot[2],
                  Mxdiff_casenew,
                  Mydiff_casenew,
                  Mzdiff_case])

    F=pd.DataFrame(F,columns=[windD])
    return F

########## Inputs Start ##########

########### Simulation input
mindist = 4
spc = 6 # maximum spacing at the beginning of the simulation for boundary creation
shape ='circle'
no_circles_squares =1 
no_turbines = 7 
file_number = f'{mindist}{mindist}{no_turbines}{no_turbines}'
windrose = 'iea'

#### Floater + turbine parameters
mass=37181396.000 #floater + turbine mass
rCG = np.array([-0.191783, 0, -11.1276])  #CG of floater + turbine
volume=36431.45 #submerged volume m3
AWP=8.514382e+02 #water plane area from .hst file 8.51441E+02*0.9996941896024465 ............8.514382e+02
rM=np.array([0, 0, -4.9540 ])

#### Mooring lines parameters 
## non permuted inputs
depth=200 #seabed depth
fairR=42.5 #fairlead radius
fair_depth=(15) #depth of fairlead

WindSpeed=10
yawang = 5
Ti = 0.06
mooringfilename = '/Turbine_' 

############################ 
path0=fr"../Results_{windrose}"
path1 =fr"../Results_{windrose}/Results_{mindist}_{spc}"
path =path1 + r"/FinalMooringWFDesign"

Aero_Forces=pd.read_csv('../Inputs/Aerodynamic_forces2.csv')

Faero=np.array([2148821.67906326,-6262.2958047442,-15538.6832881499])
Maero=np.array([17214575.6193484,3892317.0989765,2791002.68346595])

Farm_layout=pd.read_csv(path+f'/Turbines_{no_turbines}/WindSpeed_{WindSpeed}/FASTFarm_layout.csv')
xlocs  =Farm_layout.x

ylocs  = Farm_layout.y


ang = [0, 0, 0,0, 0, 0, 0, 0, 0] 
aero_angle=np.arange(0,360,5)
mydpi=96
mass=37181396.000 #floater + turbine mass
volume=36431.45 #submerged volume m3
AWP=8.514382e+02 #water plane area from .hst file 8.51441E+02*0.9996941896024465 ............8.514382e+02
rM=np.array([0, 0, -4.9540 ])


############################ steady state 
fig = plt.figure(figsize=(1000/mydpi,500/mydpi),dpi=mydpi)
# ax0 = fig.add_subplot(121)
ax1 = fig.add_subplot(111)
for i in range(0, len(xlocs)):
    

    # ms = mp.System(file=os.path.join(path,"Activefloat_MoorDyn_BS_moordyn2.dat"))
    # ms.transform(trans = [xlocs[i],ylocs[i]], rot = ang[i])
    # ms.initialize()

    # ms.plot2d(Yuvec=[0,1,0], ax=ax0,linewidth=1.5)
    # ax0.text(xlocs[i]+100,ylocs[i],f'T {i+1}',fontsize=10)
    
    ms = mp.System(file=os.path.join(path,f'Turbines_{no_turbines}',f"Turbine_{i}.dat"))
    
    ms.transform(trans = [xlocs[i],ylocs[i]], rot = ang[i])
    ms.initialize()
    100*np.cos(np.arctan(ms.lineList[0].rB[1]/ms.lineList[0].rB[0]))
    ms.plot2d(Yuvec=[0,1,0], ax=ax1,linewidth=1.5)
    
    if i==0:
        ax1.text(xlocs[i]-200,
                 ylocs[i]+200,
                 f'T {i+1}',fontsize=10)
    elif i==5:
        ax1.text(xlocs[i]-100,
                 ylocs[i]+200,
                 f'T {i+1}',fontsize=10)
    elif i==6:
        ax1.text(xlocs[i]-100,
                 ylocs[i]+200,
                 f'T {i+1}',fontsize=10)
    else:
        ax1.text(xlocs[i]+100,
                 ylocs[i]+100,
                 f'T {i+1}',fontsize=10)

ax1.set_yticklabels([])
ax1.set_xlim(-2640,2640)
# ax0.set_xlim(-2640,2640)
ax1.set_ylim(-2640,2640)
# ax0.set_ylim(-2640,2640)
ax1.set_xticks(np.arange(-2400,2880,480))
# ax0.set_xticks(np.arange(-2400,2880,480))
# ax0.set_yticks(np.arange(-2400,2880,480))
ax1.set_yticks(np.arange(-2400,2880,480))
# ax0.set_xticklabels(np.arange(-2400/240,2880/240,2))
ax1.set_xticklabels(np.arange(-2400/240,2880/240,2))
# ax0.set_yticklabels(np.arange(-2400/240,2880/240,2))
# ax0.grid()
ax1.grid()
# ax0.set_xlabel('x/D', fontdict = {'fontsize' : 12})
# ax0.set_ylabel('y/D', fontdict = {'fontsize' : 12})
# ax0.yaxis.set_label_coords(-0.1, 0.5)
ax1.set_xlabel('x/D', fontdict = {'fontsize' : 12})
# ax0.set_title('OWFL with baseline MS')
ax1.set_title('OWFL with customised MS')
# plt.subplots_adjust(wspace=0.1)
plt.show()
# pdb.set_trace()
############################
Wind_Linedict={}
for j in range(len(aero_angle)): 
    print('***************')
    print('Wind Direction ')
    print(aero_angle[j])
    print('***************')
    
    CG=np.array([-0.191783*np.cos(np.deg2rad(aero_angle[j])), -0.191783*np.sin(np.deg2rad(aero_angle[j])), -11.1276])
    # loop over mooring systems for calculating X and Xdfarm_9T_Mann_baselineiff
    Linedict = {}
    for i in range(len(xlocs)):
        # Apply the Aerodynamic forces at 8m/s and recalculate equilibrium and repeating again
        ms = mp.System(os.path.join(path,f'Turbines_{no_turbines}',f"Turbine_{i}.dat"))
        # ms = mp.System(file=os.path.join(path,"Activefloat_MoorDyn_BS_moordyn2.dat"))
        ms.transform(trans = [xlocs[i],ylocs[i]], rot = ang[i])
        
        ms.bodyList[0].rCG=CG
        ms.bodyList[0].mass=mass
        ms.bodyList[0].volume=volume
        ms.bodyList[0].AWP=AWP
        ms.bodyList[0].rM=rM
        x=np.full((6),np.nan,dtype=float)
        ms.initialize()
        # pdb.set_trace()
        x=np.full((6),np.nan,dtype=float)
        for numb in range(10):  
            if numb==0:
                F=FAero_windD(Faero,Maero,aero_angle[j],rollangle=float(0),pitchangle=float(0),yawangle=float(0))
            else:
                F=FAero_windD(Faero,Maero,aero_angle[j],rollangle=np.rad2deg(x[3]),pitchangle=np.rad2deg(x[4]),yawangle=np.rad2deg(x[5]))
            
            Forces=F[aero_angle[j]]
            ms.bodyList[0].f6Ext=Forces
            
            
            ms.solveEquilibrium(no_fail=True,maxIter=500)
            x=ms.getPositions(DOFtype="free")
        
        ms.initialize()
        if i == 0:
            fig,ax = ms.plot2d( Yuvec=[0,1,0])
          
        else:
            ms.plot2d(Yuvec=[0,1,0], ax=ax)

        L={}
        for l in range(len(ms.lineList)):
            X, Y, Z, T = ms.lineList[l].getLineCoords(1000000,n=200)
            line=LineString(np.array([X, Y, Z]).T)
            L[f'L_{l}']={'X':pd.DataFrame(X),'Y':pd.DataFrame(Y),'Z':pd.DataFrame(Z)}
        Linedict[f'T_{i}'] = L
    Wind_Linedict[f'Wind_{aero_angle[j]}']=Linedict
plt.show()


pdb.set_trace()

counter=0
summary=pd.DataFrame([], columns=['Wind_dir','Turbine_a', 'Turbine_b','min_dist'])
t_perm= list(itertools.combinations(range(no_turbines),2))
for j in range(len(aero_angle)):
    for t1 in range(len(t_perm)):
        t=t_perm[t1]
        minimum_dist=pd.DataFrame()
        for l1 in range(3):
            for l2 in range(3):
                xx=Wind_Linedict[f'Wind_{aero_angle[j]}'][f'T_{t[0]}'][f'L_{l1}']['X'].values[:, None]-Wind_Linedict[f'Wind_{aero_angle[j]}'][f'T_{t[1]}'][f'L_{l2}']['X'].values
                yy=Wind_Linedict[f'Wind_{aero_angle[j]}'][f'T_{t[0]}'][f'L_{l1}']['Y'].values[:, None]-Wind_Linedict[f'Wind_{aero_angle[j]}'][f'T_{t[1]}'][f'L_{l2}']['Y'].values
                zz=Wind_Linedict[f'Wind_{aero_angle[j]}'][f'T_{t[0]}'][f'L_{l1}']['Z'].values[:, None]-Wind_Linedict[f'Wind_{aero_angle[j]}'][f'T_{t[1]}'][f'L_{l2}']['Z'].values
                xx=pd.DataFrame(xx.reshape(len(xx),len(xx))**2)
                yy=pd.DataFrame(yy.reshape(len(xx),len(xx))**2)
                zz=pd.DataFrame(zz.reshape(len(xx),len(xx))**2)
                dist=np.sqrt(xx+yy+zz)


                summary.loc[counter,'Wind_dir'] = aero_angle[j]
                summary.loc[counter,'Turbine_a'] = f'T_{t[0]}'
                summary.loc[counter,'Turbine_b'] = f'T_{t[1]}'
                summary.loc[counter,'min_dist'] = min(dist.min())
                counter=counter+1
                
wht=summary.groupby(['Turbine_a','Turbine_b']).size()