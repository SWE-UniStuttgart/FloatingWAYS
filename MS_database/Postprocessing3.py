# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 15:34:50 2022

@author: mahfouz

This function reads the outputs from PermMatrix_fnc, cleans them up and postprocess them.
Cleaning up means it makes sure only the designs that pass all wind directions are used for calculating the displacement.
It also gives us a chance to add new constraints on the accepted designs.

We also clean up the forces, currently the forces are not used at all in any other codes. Would be nice to create a seperate function
for the forces.
We currently don't use the forces in any of the following parts of the code, but we should keep the forces post processing.

The output is saved for surge and sway motions.

The function also creates plots of the accepted mooring designs.
"""

import moorpy as mp
import numpy as np
import pandas as pd
import pdb
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.lines import Line2D
from matplotlib.animation import FuncAnimation
import random
# import os
# 

"""
I have to remove the forces and only keep the parameters I need
I can create another function only for the forces and this function only displacements
What DOFs of the floater should I keep? This is for layout optimization so maybe remove all other DOFs?
Create a function for plotting the data seperately.
Plot Headings, diameters, lengths, anchors, forces?
"""
# df=pd.DataFrame()
Xdiff_Surge=pd.DataFrame()
Xdiff_Sway=pd.DataFrame()
# Xdiff_Heave=pd.DataFrame()
# Xdiff_Roll=pd.DataFrame()
# Xdiff_Pitch=pd.DataFrame()
# Xdiff_Yaw=pd.DataFrame()
X0_Surge=pd.DataFrame()
X0_Sway=pd.DataFrame()
# X0_Yaw=pd.DataFrame()
# X0_Roll=pd.DataFrame()
X_Surge=pd.DataFrame()
X_Sway=pd.DataFrame()
# X_Heave=pd.DataFrame()
# X_Roll=pd.DataFrame()
# X_Pitch=pd.DataFrame()
# X_Yaw=pd.DataFrame()
F0X1=pd.DataFrame()
F0X2=pd.DataFrame()
F0X3=pd.DataFrame()
FX1=pd.DataFrame()
FX2=pd.DataFrame()
FX3=pd.DataFrame()
F0Y1=pd.DataFrame()
F0Y2=pd.DataFrame()
F0Y3=pd.DataFrame()
FY1=pd.DataFrame()
FY2=pd.DataFrame()
FY3=pd.DataFrame()
# F_inline1= pd.DataFrame()
# F_perp1= pd.DataFrame()
# F_inline2=pd.DataFrame()
# F_perp2=pd.DataFrame()
# F_inline3=pd.DataFrame()
# F_perp3= pd.DataFrame()

path= r"D:\Moordyn\moorpy\MoorPy\Outputs\trynew\WindSpeed_10"
path0= r"D:\Moordyn\moorpy\MoorPy\Outputs\trynew"


########### Reading the outputs only taking the accepted designs 
LastPermMatrix= pd.read_csv(path+ '\PermutationMatrix355.csv',header=0, names=['PermNum','ang1','D1','ancR1','L1','ang2','D2','ancR2','L2','ang3','D3','ancR3','L3'])
# LastPermMatrix= pd.read_csv(path+ '\PermutationMatrix350.csv',header=0, names=['PermNum','D1', 'ang1', 'ancR1','L1','D2', 'ang2', 'ancR2','L2','D3', 'ang3', 'ancR3','L3'])
AcceptedPerm=pd.DataFrame(LastPermMatrix.PermNum)
aero_angle=np.arange(0,360,5)
for i in range(len(aero_angle)): 
    if i==0:
        PermMatrix= pd.read_csv(path0+'\PermutationMatrix.csv',header=0, names=['PermNum','ang1','D1','ancR1','L1','ang2','D2','ancR2','L2','ang3','D3','ancR3','L3'])
        # PermMatrix= pd.read_csv(path+'\PermutationMatrix.csv',header=0, names=['PermNum','D1', 'ang1', 'ancR1','L1','D2', 'ang2', 'ancR2','L2','D3', 'ang3', 'ancR3','L3'])
    else:
        PermMatrix= pd.read_csv(path+'\PermutationMatrix'+str(aero_angle[i-1])+'.csv',header=0, names=['PermNum','ang1','D1','ancR1','L1','ang2','D2','ancR2','L2','ang3','D3','ancR3','L3'])
        # PermMatrix= pd.read_csv(path+'\PermutationMatrix'+str(aero_angle[i-1])+'.csv',header=0, names=['PermNum','D1', 'ang1', 'ancR1','L1','D2', 'ang2', 'ancR2','L2','D3', 'ang3', 'ancR3','L3'])
    PermMatrix=PermMatrix.loc[PermMatrix['PermNum'].isin(AcceptedPerm.PermNum)]
    Results = pd.read_csv(path+'\FinalPermutationResult'+str(aero_angle[i])+'.csv',header=0,names=['displacement','PermNum','Surge','Sway','Heave','Roll','Pitch','Yaw'])
    # pdb.set_trace()
    Results= Results.loc[Results['PermNum'].isin(AcceptedPerm.PermNum)]
    # pdb.set_trace()
    Forces = pd.read_csv(path+'\ForceResult'+str(aero_angle[i])+'.csv',header=0,names=['direction','PermNum','F0L1','F0L2','F0L3','FL1','FL2','FL3'])
    Forces = Forces.loc[Forces['PermNum'].isin(AcceptedPerm.PermNum)]
    ForceX = pd.DataFrame(Forces[(Forces['direction'] == 'X')])
    ForceX.reset_index(drop=True, inplace=True)
    F0X1[aero_angle[i]]=ForceX.F0L1.values
    F0X2[aero_angle[i]]=ForceX.F0L2.values
    F0X3[aero_angle[i]]=ForceX.F0L3.values
    FX1[aero_angle[i]]=ForceX.FL1.values
    FX2[aero_angle[i]]=ForceX.FL2.values
    FX3[aero_angle[i]]=ForceX.FL3.values
    ForceY = pd.DataFrame(Forces[(Forces['direction'] == 'Y')])
    ForceY.reset_index(drop=True, inplace=True)
    F0Y1[aero_angle[i]]=ForceY.F0L1.values
    F0Y2[aero_angle[i]]=ForceY.F0L2.values
    F0Y3[aero_angle[i]]=ForceY.F0L3.values
    FY1[aero_angle[i]]=ForceY.FL1.values
    FY2[aero_angle[i]]=ForceY.FL2.values
    FY3[aero_angle[i]]=ForceY.FL3.values
    Xdiff=pd.DataFrame(Results[(Results['displacement'] == 'X_diff')])
    Xdiff.reset_index(drop=True, inplace=True)
    Xdiff_Surge[aero_angle[i]]=Xdiff.Surge.values
    Xdiff_Sway[aero_angle[i]]=Xdiff.Sway.values
    # Xdiff_Heave[aero_angle[i]]=Xdiff.Heave.values
    # Xdiff_Roll[aero_angle[i]]=Xdiff.Roll.values*180/np.pi
    # Xdiff_Pitch[aero_angle[i]]=Xdiff.Pitch.values*180/np.pi
    # Xdiff_Yaw[aero_angle[i]]=Xdiff.Yaw.values*180/np.pi
    X=pd.DataFrame(Results[(Results['displacement'] == 'X')])
    X.reset_index(drop=True, inplace=True)
    X_Surge[aero_angle[i]]=X.Surge.values
    X_Sway[aero_angle[i]]=X.Sway.values
    # X_Heave[aero_angle[i]]=X.Heave.values
    # X_Roll[aero_angle[i]]=X.Roll.values*180/np.pi
    # X_Pitch[aero_angle[i]]=X.Pitch.values*180/np.pi
    # X_Yaw[aero_angle[i]]=X.Yaw.values*180/np.pi
    X0=pd.DataFrame(Results[(Results['displacement'] == 'X0')])
    X0.reset_index(drop=True, inplace=True)
    X0_Surge[aero_angle[i]]=X0.Surge.values
    X0_Sway[aero_angle[i]]=X0.Sway.values
    # X0_Yaw[aero_angle[i]]=X0.Yaw.values*180/np.pi
    # X0_Roll[aero_angle[i]]=X0.Roll.values*180/np.pi
    PermMatrix.reset_index(drop=True, inplace=True)
########### end of reading the outputs only taking the accepted designs    


########### Rotating the refernce frame. 
########### To make surge always in wind direction motion and sway perpendicular to it
########### This is repeated in other functions so should be a function on its own.
df=pd.DataFrame(np.tile(-aero_angle,(len(PermMatrix),1)), columns=np.arange(0,360,5))
Xnew_Surge=X_Surge*np.cos(np.deg2rad(-df))+X_Sway*np.sin(np.deg2rad(-df))
Xnew_Sway=X_Sway*np.cos(np.deg2rad(-df))-X_Surge*np.sin(np.deg2rad(-df))

X0new_Surge=X0_Surge*np.cos(np.deg2rad(-df))+X0_Sway*np.sin(np.deg2rad(-df))
X0new_Sway=X0_Sway*np.cos(np.deg2rad(-df))-X0_Surge*np.sin(np.deg2rad(-df))

########### end otating the refernce frame. 

Xdiffnew_Surge=np.absolute(Xnew_Surge-X0new_Surge)
Xdiffnew_Sway=np.absolute(Xnew_Sway-X0new_Sway)
# Xdiffnew_Surge=Xnew_Surge-X0new_Surge
# Xdiffnew_Sway=Xnew_Sway-X0new_Sway

Logic=np.zeros(len(Xdiffnew_Surge))
for j in range(len(Xdiffnew_Surge)):
    Logic[j]=np.max(Xdiffnew_Surge.iloc[j,:]) <= 240
    
# Logic=np.zeros(len(Xdiffnew_Sway))
# for j in range(len(Xdiffnew_Sway)):
#     Logic[j]=np.max(Xdiffnew_Sway.iloc[j,:]) >= 0


########### Adding new constraint
########### This is a bad way to do it I should change this
Logic=pd.DataFrame(Logic)
MinimumPerpDisp=pd.DataFrame(Logic[(Logic[0] == 1)])

FinalPermMatrix=PermMatrix.loc[PermMatrix.index.isin(MinimumPerpDisp.index)]
FinalPermMatrix.reset_index(drop=True, inplace=True)

Xnew_Sway=Xnew_Sway.loc[Xnew_Sway.index.isin(MinimumPerpDisp.index)]  
Xnew_Sway.reset_index(drop=True, inplace=True)

Xnew_Surge=Xnew_Surge.loc[Xnew_Surge.index.isin(MinimumPerpDisp.index)]  
Xnew_Surge.reset_index(drop=True, inplace=True)

X0new_Surge=X0new_Surge.loc[X0new_Surge.index.isin(MinimumPerpDisp.index)]  
X0new_Surge.reset_index(drop=True, inplace=True)

X0new_Sway=X0new_Sway.loc[X0new_Sway.index.isin(MinimumPerpDisp.index)]  
X0new_Sway.reset_index(drop=True, inplace=True)

Xdiff_Surge=Xdiff_Surge.loc[Xdiff_Surge.index.isin(MinimumPerpDisp.index)]  
Xdiff_Surge.reset_index(drop=True, inplace=True)

Xdiff_Sway=Xdiff_Sway.loc[Xdiff_Sway.index.isin(MinimumPerpDisp.index)]  
Xdiff_Sway.reset_index(drop=True, inplace=True)

Xdiffnew_Surge=Xdiffnew_Surge.loc[Xdiffnew_Surge.index.isin(MinimumPerpDisp.index)]  
Xdiffnew_Surge.reset_index(drop=True, inplace=True)

Xdiffnew_Sway=Xdiffnew_Sway.loc[Xdiffnew_Sway.index.isin(MinimumPerpDisp.index)]  
Xdiffnew_Sway.reset_index(drop=True, inplace=True)

X0_Surge=X0_Surge.loc[X0_Surge.index.isin(MinimumPerpDisp.index)]  
X0_Surge.reset_index(drop=True, inplace=True)

X0_Sway=X0_Sway.loc[X0_Sway.index.isin(MinimumPerpDisp.index)]  
X0_Sway.reset_index(drop=True, inplace=True)

X_Surge=X_Surge.loc[X_Surge.index.isin(MinimumPerpDisp.index)]  
X_Surge.reset_index(drop=True, inplace=True)

X_Sway=X_Sway.loc[X_Sway.index.isin(MinimumPerpDisp.index)] 
X_Sway.reset_index(drop=True, inplace=True)

F0X1=F0X1.loc[F0X1.index.isin(MinimumPerpDisp.index)] 
F0X1.reset_index(drop=True, inplace=True)

F0X2=F0X2.loc[F0X2.index.isin(MinimumPerpDisp.index)] 
F0X2.reset_index(drop=True, inplace=True)

F0X3=F0X3.loc[F0X3.index.isin(MinimumPerpDisp.index)] 
F0X3.reset_index(drop=True, inplace=True)

FX1=FX1.loc[FX1.index.isin(MinimumPerpDisp.index)] 
FX1.reset_index(drop=True, inplace=True) 

FX2=FX2.loc[FX2.index.isin(MinimumPerpDisp.index)] 
FX2.reset_index(drop=True, inplace=True) 

FX3=FX3.loc[FX3.index.isin(MinimumPerpDisp.index)] 
FX3.reset_index(drop=True, inplace=True) 

F0Y1=F0Y1.loc[F0Y1.index.isin(MinimumPerpDisp.index)] 
F0Y1.reset_index(drop=True, inplace=True) 

F0Y2=F0Y2.loc[F0Y2.index.isin(MinimumPerpDisp.index)] 
F0Y2.reset_index(drop=True, inplace=True) 

F0Y3=F0Y3.loc[F0Y3.index.isin(MinimumPerpDisp.index)] 
F0Y3.reset_index(drop=True, inplace=True) 

FY1=FY1.loc[FY1.index.isin(MinimumPerpDisp.index)] 
FY1.reset_index(drop=True, inplace=True) 

FY2=FY2.loc[FY2.index.isin(MinimumPerpDisp.index)] 
FY2.reset_index(drop=True, inplace=True) 

FY3=FY3.loc[FY3.index.isin(MinimumPerpDisp.index)] 
FY3.reset_index(drop=True, inplace=True) 
########### end adding new constraint

FX_Sum=pd.DataFrame(np.abs(FX1+FX2+FX3))
FY_Sum=pd.DataFrame(np.abs(FY1+FY2+FY3))
F0X_Sum=pd.DataFrame(np.abs(F0X1+F0X2+F0X3))
F0Y_Sum=pd.DataFrame(np.abs(F0Y1+F0Y2+F0Y3))
F_resultant=pd.DataFrame(np.sqrt(FX_Sum**2+FY_Sum**2))
F0_resultant=pd.DataFrame(np.sqrt(F0X_Sum**2+F0Y_Sum**2))
F1Res=pd.DataFrame(np.sqrt(FX1**2+FY1**2))
F2Res=pd.DataFrame(np.sqrt(FX2**2+FY2**2))
F3Res=pd.DataFrame(np.sqrt(FX3**2+FY3**2))
F0Res1=pd.DataFrame(np.sqrt(F0X1**2+F0Y1**2))
F0Res2=pd.DataFrame(np.sqrt(F0X2**2+F0Y2**2))
F0Res3=pd.DataFrame(np.sqrt(F0X3**2+F0Y3**2))
F1ResAngle=np.rad2deg(np.arctan2(FY1, FX1))
F2ResAngle=np.rad2deg(np.arctan2(FY2, FX2))
F3ResAngle=np.rad2deg(np.arctan2(FY3, FX3))+360
F01ResAngle=np.rad2deg(np.arctan2(F0Y1, F0X1))
F02ResAngle=np.rad2deg(np.arctan2(F0Y2, F0X2))
F03ResAngle=np.rad2deg(np.arctan2(F0Y3, F0X3))+360


df=pd.DataFrame(np.tile(-aero_angle,(len(FinalPermMatrix),1)), columns=np.arange(0,360,5))
F_inline1= F1Res.values*np.cos(np.deg2rad(df.sub(-F1ResAngle,axis='index')))
F_perp1=F1Res.values*np.sin(np.deg2rad(df.sub(-F1ResAngle,axis='index')))
F_inline2= F2Res.values*np.cos(np.deg2rad(df.sub(-F2ResAngle,axis='index')))
F_perp2= F2Res.values*np.sin(np.deg2rad(df.sub(-F2ResAngle,axis='index')))
F_inline3= F3Res.values*np.cos(np.deg2rad(df.sub(-F3ResAngle,axis='index')))
F_perp3= F3Res.values*np.sin(np.deg2rad(df.sub(-F3ResAngle,axis='index')))


F0_inline1= F0Res1.values*np.cos(np.deg2rad(df.sub(-F01ResAngle,axis='index')))
F0_perp1=F0Res1.values*np.sin(np.deg2rad(df.sub(-F01ResAngle,axis='index')))
F0_inline2= F0Res2.values*np.cos(np.deg2rad(df.sub(-F02ResAngle,axis='index')))
F0_perp2= F0Res2.values*np.sin(np.deg2rad(df.sub(-F02ResAngle,axis='index')))
F0_inline3= F0Res3.values*np.cos(np.deg2rad(df.sub(-F03ResAngle,axis='index')))
F0_perp3= F0Res3.values*np.sin(np.deg2rad(df.sub(-F03ResAngle,axis='index')))




############  start: Calculating the difference in displacement assuming no surge motion only sway
#  Getting the x and y coordinates assuming only sway motion
deltaX_perp=X_Surge- (X0_Surge+Xdiffnew_Surge*np.cos(np.deg2rad(-df)))
deltaY_perp=X_Sway-(X0_Sway+Xdiffnew_Surge*np.sin(np.deg2rad(-df)))
############ end ############

perp_disp=Xnew_Sway-X0new_Sway

# deltaX_perp2=X_Surge- (X0_Surge+(Xnew_Surge-X0new_Surge)*np.cos(np.deg2rad(-df)))
# deltaY_perp2=X_Sway-(X0_Sway+(Xnew_Surge-X0new_Surge)*np.sin(np.deg2rad(-df)))
# Xdiffnew_Sway=Xnew_Sway-X0new_Sway

############ These two files are not used afterwards ############
deltaX_perp.to_csv(path+ '\Position_X_noinline'+'.csv', index=False)
deltaY_perp.to_csv(path+ '\Position_Y_noinline'+'.csv', index=False)
############ end ############

FinalPermMatrix.to_csv(path+ '\FinalPermMatrix'+'.csv', index=False)
Xdiff_Surge.to_csv(path+ '\Xdiff'+'.csv', index=False)
Xdiff_Sway.to_csv(path+ '\Ydiff'+'.csv', index=False)
perp_disp.to_csv(path+ '\perpdiff'+'.csv', index=False)

pdb.set_trace() 

########## The plotiing part starts

Angles=pd.DataFrame(FinalPermMatrix,columns=['ang1','ang2','ang3'])
Angles.drop_duplicates(inplace=True, ignore_index=True)
Diameters=pd.DataFrame(FinalPermMatrix,columns=['D1','D2','D3'])
Diameters.drop_duplicates(inplace=True, ignore_index=True)
Length=pd.DataFrame(FinalPermMatrix,columns=['L1','ancR1','L2','ancR2','L3','ancR3'])
Length.drop_duplicates(inplace=True, ignore_index=True)
AncR=pd.DataFrame(FinalPermMatrix,columns=['ancR1','ancR2','ancR3'])
AncR.drop_duplicates(inplace=True, ignore_index=True)

Perm_Angles = {}
Perm_Angles_Diameters = {}
Perm_Angles_Length = {}
Perm_Angles_ancR = {}
for n in range(len(Angles)):
    Perm_Angles['Angles_%s_%s_%s' %(Angles['ang1'][n],  Angles['ang2'][n], Angles['ang3'][n])] = pd.DataFrame(FinalPermMatrix[(FinalPermMatrix['ang1'] == Angles['ang1'][n]) & (FinalPermMatrix['ang2'] ==Angles['ang2'][n]) & (FinalPermMatrix['ang3'] == Angles['ang3'][n])])

key=list(Perm_Angles.keys())
for n in range(len(Perm_Angles.keys())):
    for i in range(len(Diameters)):
        Perm_Angles_Diameters[ key[n] +'_Diameters_%s_%s_%s' %(Diameters['D1'][i],  Diameters['D2'][i], Diameters['D3'][i])] = Perm_Angles[key[n]].loc[(Perm_Angles[key[n]].D1==Diameters['D1'][i] ) & (Perm_Angles[key[n]].D2==Diameters['D2'][i] ) & (Perm_Angles[key[n]].D3== Diameters['D3'][i])]
# Perm_Angles_Diameters={k:v for (k,v) in Perm_Angles_Diameters.items() if not v.empty}        
for n in range(len(Perm_Angles.keys())):
    for i in range(len(Length)):
        Perm_Angles_Length[ key[n] +'_Length_%s_%s_%s' %(Length['L1'][i],  Length['L2'][i], Length['L3'][i])] = Perm_Angles[key[n]].loc[(Perm_Angles[key[n]].L1==Length['L1'][i] ) & (Perm_Angles[key[n]].L2==Length['L2'][i] ) & (Perm_Angles[key[n]].L3== Length['L3'][i])]
# Perm_Angles_Length={k:v for (k,v) in Perm_Angles_Length.items() if not v.empty}
for n in range(len(Perm_Angles.keys())):
    for i in range(len(AncR)):
        Perm_Angles_ancR[ key[n] +'_AncR_%s_%s_%s' %(AncR['ancR1'][i],  AncR['ancR2'][i], AncR['ancR3'][i])] = Perm_Angles[key[n]].loc[(Perm_Angles[key[n]].ancR1==AncR['ancR1'][i] ) & (Perm_Angles[key[n]].ancR2==AncR['ancR2'][i] ) & (Perm_Angles[key[n]].ancR3== AncR['ancR3'][i])]
# Perm_Angles_ancR={k:v for (k,v) in Perm_Angles_ancR.items() if not v.empty}


pdb.set_trace() 
   

#  paper figures

def add_scale(ax):
    # add extra axes for the scale
    rect = ax.get_position()
    rect = (rect.xmin+X_OFFSET, rect.ymin+rect.height/2+Y_OFFSET, # x, y
            rect.width, rect.height/2) # width, height
    scale_ax = ax.figure.add_axes(rect)
    # hide most elements of the new axes
    for loc in ['right', 'top', 'bottom']:
        scale_ax.spines[loc].set_visible(False)
    scale_ax.tick_params(bottom=False, labelbottom=False)
    scale_ax.patch.set_visible(False) # hide white background
    # adjust the scale
    scale_ax.spines['left'].set_bounds(*ax.get_ylim())
    # scale_ax.spines['left'].set_bounds(0, ax.get_rmax()) # mpl < 2.2.3
    scale_ax.set_yticks(ax.get_yticks())
    scale_ax.set_ylim(ax.get_rorigin(), ax.get_rmax())
    scale_ax.set_yticklabels(ax.get_yticks(), fontdict = {'fontsize' : 14})
    # scale_ax.set_ylim(ax.get_ylim()) # Matplotlib < 2.2.3
for m in range(len(Angles)):
# for m in range(0,28,4):
    # m=15
    Perm= Perm_Angles['Angles_%s_%s_%s' %(Angles['ang1'][m],  Angles['ang2'][m], Angles['ang3'][m])]
    colors=['r','b','k','g','y']
    mydpi=96
    fig = plt.figure(figsize=(1500/mydpi,550/mydpi),dpi=mydpi)
    ax0 = fig.add_subplot(131)
    ax1 = fig.add_subplot(132,projection='polar')
    ax2 = fig.add_subplot(133,projection='polar')
    counter=0
    for i in random.sample(range(0, len(Perm)), 1):
        
        ax0.plot(np.append(Xdiff_Surge.iloc[Perm.index[i]],Xdiff_Surge.iloc[Perm.index[i],0]), 
                 np.append(Xdiff_Sway.iloc[Perm.index[i]],Xdiff_Sway.iloc[Perm.index[i],0]),colors[counter]+'.')
        ax1.plot(np.append(np.arange(0,2.0*np.pi,np.pi/36),0), 
                 np.append(Xdiffnew_Sway.iloc[Perm.index[i]].values, np.array(Xdiffnew_Sway.iloc[Perm.index[i],0])),colors[counter])
        # ax.plot(np.append(np.arange(0,2.0*np.pi,np.pi/18.0),0), np.append(Angle12.iloc[Perm.index[i]].values, np.array(Angle12.iloc[Perm.index[i],0])),'-g')
        ax1.arrow(np.deg2rad(Perm['ang1'].iloc[i]), 0, 0, 240, head_width=0, head_length=0, alpha = 1, width = 0.0,
                  edgecolor = 'black', facecolor = 'black', zorder = 2)
        ax1.arrow(np.deg2rad(Perm['ang2'].iloc[i]), 0, 0, 240, head_width=0, head_length=0, alpha = 1, width = 0.0,
                  edgecolor = 'black', facecolor = 'black',  zorder = 2)
        ax1.arrow(np.deg2rad(Perm['ang3'].iloc[i]), 0, 0, 240, head_width=0, head_length=0, alpha = 1, width = 0.0,
                  edgecolor = 'black', facecolor = 'black',  zorder = 2)
        ax2.plot(np.append(np.arange(0,2.0*np.pi,np.pi/36),0), 
                 np.append(Xdiffnew_Surge.iloc[Perm.index[i]].values, np.array(Xdiffnew_Surge.iloc[Perm.index[i],0])),colors[counter])
        ax2.arrow(np.deg2rad(Perm['ang1'].iloc[i]), 0, 0, 240, head_width=0, head_length=0, alpha = 1, width = 0.0,
                  edgecolor = 'black', facecolor = 'black',  zorder = 2)
        ax2.arrow(np.deg2rad(Perm['ang2'].iloc[i]), 0, 0, 240, head_width=0, head_length=0, alpha = 1, width = 0.0,
                  edgecolor = 'black', facecolor = 'black',  zorder = 2)
        ax2.arrow(np.deg2rad(Perm['ang3'].iloc[i]), 0, 0, 240, head_width=0, head_length=0, alpha = 1, width = 0.0,
                  edgecolor = 'black', facecolor = 'black', zorder = 2)
        ax0.set_xticks(np.arange(-300,400,100))
        ax0.set_yticks(np.arange(-300,400,100))
        ax0.set_xticklabels(np.arange(-300,400,100), fontdict = {'fontsize' : 14})
        ax0.set_yticklabels(np.arange(-300,400,100), fontdict = {'fontsize' : 14})
        # ax1.set_rticks([50, 100, 150, 200,250,300,350])  # Less radial ticks
        ax1.set_xticks(np.arange(0,2.0*np.pi,np.pi/9.0))
        ax1.set_xticklabels(np.arange(0,360,20), fontdict = {'fontsize' : 14})
        ax1.set_yticklabels([])
        ax1.set_rticks([50,100,150,200,250])  # Less radial ticks
        ax1.set_xticks(np.arange(0,2.0*np.pi,np.pi/9.0))
        ax1.set_xticklabels(np.arange(0,360,20), fontdict = {'fontsize' : 14})
        ax1.set_yticklabels([])
        ax2.set_rticks([50,100,150,200,250])  # Less radial ticks
        ax2.set_xticks(np.arange(0,2.0*np.pi,np.pi/9.0))
        ax2.set_xticklabels(np.arange(0,360,20), fontdict = {'fontsize' : 14})
        ax2.set_yticklabels([])

        counter=counter+1
    ax0.set_title('FOWT watch circle \n Angles_%s_%s_%s' %(Angles['ang1'][m],  Angles['ang2'][m], Angles['ang3'][m]), fontdict = {'fontsize' : 18})
    ax0.grid(True)
    ax0.set_xlabel('X [m]', fontdict = {'fontsize' : 14})
    ax0.set_ylabel('Y [m]', fontdict = {'fontsize' : 14})
    X_OFFSET = -0.02
    Y_OFFSET = 0.11
    add_scale(ax2)
    ax1.set_title('Perp displacement[m] \n Angles_%s_%s_%s' %(Angles['ang1'][m],  Angles['ang2'][m], Angles['ang3'][m]), fontdict = {'fontsize' : 18})
    ax2.set_title('Inline displacement [m] \n Angles_%s_%s_%s' %(Angles['ang1'][m],  Angles['ang2'][m], Angles['ang3'][m]), fontdict = {'fontsize' : 18})
    # plt.savefig(r'D:\data\40_mahfouz\Moordyn\MoorPy\Torque2022\Figures'+ '\Angles_%s_%s_%s.png' %(Angles['ang1'][m],  Angles['ang2'][m], Angles['ang3'][m]),dpi=200)
    plt.show()