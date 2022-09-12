# -*- coding: utf-8 -*-
"""
Created on Thu Sep  9 17:23:43 2021


This function creates the permutation matrix and then checks the mooring lines that pass the defined constraints.

The constraints of all degrees of freedom except the constraints on surge and sway motion can be done here.

The surge and swat motion constraints can only be done in the post processing.

We should make the constraint about the vertical forces on the anchor optionsl.

Step1: Creating the permutation matrix
Step2: loop over wind directions
Step3: loop over mooring designs in the permutation matrix
Step4: Calculating static equilibrium
Step5: Applying the aerodynamic force
Step6: Calculating the new equilibrium
@author: mahfouz
"""

import numpy as np
from itertools import product 
import pandas as pd
import moorpy as mp
import pdb
import os



def PermutationMatrix(depth, fairR, D, angle1, angle2, angle3, ancr, L_ratio):
    D=np.array(D)
    angle1=np.array(angle1)
    angle2=np.array(angle2)
    angle3=np.array(angle3)
    ancr=np.array(ancr)
    ancA_L=[]



    permute=[angle1,angle2,angle3]
    Angles=pd.DataFrame(list(product(*permute)), columns=['ang1','ang2', 'ang3'])
    Angles_diff=pd.DataFrame()
    Angles_diff['Ang21_diff'] = Angles['ang2'] - Angles['ang1']
    Angles_diff['Ang32_diff'] = Angles['ang3'] - Angles['ang2']
    Angles_diff['Ang31_diff'] = np.absolute(Angles['ang3'] - Angles['ang1']-360)
    rpeated= Angles_diff[pd.DataFrame(np.sort(Angles_diff.values), columns=Angles_diff.columns, index=Angles_diff.index).duplicated(keep='first')]
    Angles=Angles.loc[~Angles.index.isin(rpeated.index)]
    Angles=Angles.values.tolist()
## This loop is to find the permutation between lines length and anchors radius 
    if np.size(ancr)==1:

        B=list(product(*[np.array([ancr]),np.array(L_ratio)])) #creates the relation between ancr and length should be done in a smarter way
        ancA_L.extend(B)    
    else:       
        for i in range(len(ancr)):
            B=list(product(*[np.array([ancr[i]]),np.array(L_ratio)])) #creates the relation between ancr and length should be done in a smarter way
            ancA_L.extend(B)

    permute=[Angles,D,ancA_L,D,ancA_L,D,ancA_L]
    # Create permutation matrix using pandas dataframes
    PermutMatrix=pd.DataFrame(list(product(*permute)), columns=['angles','D1', 'anc1','D2',  'anc2','D3',  'anc3']) 
    # Lmin=np.round_(np.sqrt((ancr-fairR)**2+depth**2),decimals=1)
    # Lmax=ancr+depth-fairR
    # df2=pd.DataFrame(PermutMatrix)
    df2=pd.DataFrame(PermutMatrix.anc1.tolist(), index= PermutMatrix.index)
    # pdb.set_trace()
    PermutMatrix['ancR1']=df2[0]
    Lmin=np.round_(np.sqrt((df2[0]-fairR)**2+depth**2),decimals=1)
    Lmax=df2[0]+depth-fairR
    PermutMatrix['L1']=Lmin+df2[1]*(Lmax-Lmin)
    df2=pd.DataFrame(PermutMatrix.anc2.tolist(), index= PermutMatrix.index)
    PermutMatrix['ancR2']=df2[0]
    Lmin=np.round_(np.sqrt((df2[0]-fairR)**2+depth**2),decimals=1)
    Lmax=df2[0]+depth-fairR
    PermutMatrix['L2']=Lmin+df2[1]*(Lmax-Lmin)
    df2=pd.DataFrame(PermutMatrix.anc3.tolist(), index= PermutMatrix.index)
    PermutMatrix['ancR3']=df2[0]
    Lmin=np.round_(np.sqrt((df2[0]-fairR)**2+depth**2),decimals=1)
    Lmax=df2[0]+depth-fairR
    PermutMatrix['L3']=Lmin+df2[1]*(Lmax-Lmin)
    df2=pd.DataFrame(PermutMatrix.angles.tolist(), index= PermutMatrix.index)
    PermutMatrix['ang1']=df2[0]
    PermutMatrix['ang2']=df2[1]
    PermutMatrix['ang3']=df2[2]
    PermutMatrix['PermNum']=PermutMatrix.index
    PermMatrix=pd.DataFrame(PermutMatrix,columns=['PermNum','ang1', 'D1', 'ancR1','L1','ang2', 'D2', 'ancR2','L2','ang3', 'D3', 'ancR3','L3'])

    del PermutMatrix
    del df2
    del B

    # making sure no lines have the same anchor angle
    PermMatrix = PermMatrix.drop(PermMatrix[PermMatrix.ang1==PermMatrix.ang2].index)
    PermMatrix = PermMatrix.drop(PermMatrix[PermMatrix.ang2==PermMatrix.ang3].index)
    PermMatrix = PermMatrix.drop(PermMatrix[PermMatrix.ang1==PermMatrix.ang3-360].index)
    PermMatrix.reset_index(drop=True, inplace=True)
    #Save the permutation matrix as csv
    # PermMatrix.to_csv('PermutationMatrix.csv')
    return PermMatrix

def FAero_windD(Faero,windD):
    #  Decomposes the aero dynamic force in x and y components according to the wind direction
    F=np.array([Faero[0]*np.cos(np.deg2rad(windD))+Faero[1]*np.sin(np.deg2rad(windD)),
                 Faero[1]*np.cos(np.deg2rad(windD))+Faero[0]*np.sin(np.deg2rad(windD)),
                 np.repeat(Faero[2],len(windD)),
                 Faero[3]*np.cos(np.deg2rad(windD))+Faero[4]*np.sin(np.deg2rad(windD)),
                 Faero[4]*np.cos(np.deg2rad(windD))+Faero[3]*np.sin(np.deg2rad(windD)),
                 np.repeat(Faero[5],len(windD))])
    F=pd.DataFrame(F,columns=windD)
    return F

def makeSampleSystem(depth, type_string, LineD, dryMass_L, EA, angle, anchorR, fair_depth, fairR, LineLength,mass,volume,rCG,AWP,rM):
    #  This function creates the mooring system from the mooring design inputs
    #depth: integer, seabed depth from seawater level positive [m]
    #type_string: string, name of each line to define it's characteristic
    #LineD: diameter of each line [m]
    #dryMass_: dry mass of the line [kg/m]
    #EA: axial stiffness [N]
    #angle: angle of the fairlead and anchor [deg]
    #anchorR: anchor radius [m]
    #fair_depth: depth of the fairlead from sea water level [m]
    #fairR: fairlead radius [m]
    #LineLength: length of each line [m]
    #mass: mass of floater
    #volume: volume of floater
    #rCG: CG if the floater
    #AWP: water plane area
    #rm: metacenter relative to the floater 
    
    # Create blank system object
    MooringSystem = mp.System()
    
    # Set the depth of the system to the depth of the input value
    MooringSystem.depth = depth
    
    # For the number of variables in the LineD variable, add each Linetype to the system
    for i in range(len(type_string)):
        MooringSystem.addLineType(type_string[i], LineD[i], dryMass_L[i], EA[i])
    
    # Add a blank, free, body at [0,0,0] to the system
    MooringSystem.addBody(0, np.zeros(6),mass,volume,rCG,AWP,rM)
    
    # Set the anchor points of the system
    anchors = []
    for i in range(len(angle)):
        MooringSystem.addPoint(1, np.array([anchorR[i]*np.cos(np.deg2rad(angle[i])), anchorR[i]*np.sin(np.deg2rad(angle[i])), -MooringSystem.depth], dtype=float))
        anchors.append(len(MooringSystem.pointList))

    # Set the points that are attached to the body to the system
    bodypts = []
    for i in range(len(angle)):
        MooringSystem.addPoint(1, np.array([fairR*np.cos(np.deg2rad(angle[i])), fairR*np.sin(np.deg2rad(angle[i])), -fair_depth], dtype=float))
        bodypts.append(len(MooringSystem.pointList))
        MooringSystem.bodyList[0].attachPoint(MooringSystem.pointList[bodypts[i]-1].number, MooringSystem.pointList[bodypts[i]-1].r - MooringSystem.bodyList[0].r6[:3])
    
    # Add and attach lines to go from the anchor points to the body points
    for i in range(len(angle)):    
        MooringSystem.addLine(LineLength[i], type_string[i], nSegs=80)
        line = len(MooringSystem.lineList)
        MooringSystem.pointList[anchors[i]-1].attachLine(MooringSystem.lineList[line-1].number, 0)
        MooringSystem.pointList[bodypts[i]-1].attachLine(MooringSystem.lineList[line-1].number, 1)
        
    return MooringSystem

def MotionConstraints (Const_dict, path, resultfile,PermMatrix):
    # This function doesnt only drop constraints in the defined DoFs.
    # It also drops the mooring system that didnt achieve the other constraints or MoorPy failed to find their equilibrium position
    Results= pd.read_csv(path+resultfile,header=0,names=['displacement','permnum','Surge','Sway','Heave','Roll','Pitch','Yaw'])
    df=Results
    dd=pd.DataFrame(df[((df['displacement'] == 'X_diff') & (df['Surge'] == 0) & (df['Sway'] == 0) & (df['Heave'] == 0) & (df['Roll'] == 0) & (df['Pitch'] == 0) & (df['Yaw'] == 0))].permnum)
    df=df.loc[~df.permnum.isin(dd.permnum)]
    dd=pd.DataFrame(df[((df['displacement'] == 'X_diff') & (df['Surge'] == 0) & (df['Sway'] == 0) & (df['Yaw'] == 0))].permnum)
    df=df.loc[~df.permnum.isin(dd.permnum)]
    df.Roll=df.Roll*180/np.pi
    df.Pitch=df.Pitch*180/np.pi
    df.Yaw=df.Yaw*180/np.pi
    inputs=pd.DataFrame(Const_dict)
    col=inputs.columns
    for i in range(len(col)):
        dd=pd.DataFrame(df[(df['displacement'] == 'X_diff') & ((df[col[i]] > inputs.iloc[2,i][1]) | (df[col[i]] < inputs.iloc[2,i][0]))].permnum)
        df=df.loc[~df.permnum.isin(dd.permnum)]
        dd=pd.DataFrame(df[(df['displacement'] == 'X0') & ((df[col[i]] > inputs.iloc[0,i][1]) | (df[col[i]] < inputs.iloc[0,i][0]))].permnum)
        df=df.loc[~df.permnum.isin(dd.permnum)]
        dd=pd.DataFrame(df[(df['displacement'] == 'X') & ((df[col[i]] > inputs.iloc[1,i][1]) | (df[col[i]] < inputs.iloc[1,i][0]))].permnum)
        df=df.loc[~df.permnum.isin(dd.permnum)]

    PermNum=pd.DataFrame(np.unique(df.permnum))
    # PermNum.reset_index(drop=True, inplace=True)
    # pdb.set_trace()
    PermMatrix=PermMatrix.loc[PermMatrix.PermNum.isin(PermNum[0])]
    PermMatrix.reset_index(drop=True, inplace=True)
    return PermMatrix


print('Creating Permutation matrix')
########## Inputs Start ##########
#### Floater + turbine parameters
mass=3.66445e+07 #floater + turbine mass
rCG=np.array([-0.194593, 0, -11.0362]) #CG of floater + turbine
volume=36431.45 #submerged volume m3
AWP=8.51441E+02 #water plane area from AQWA
rM=np.array([0, 0, -4.9540]) #metacentric height from AQWA (CG - metacentricheight)

#### Mooring lines parameters 
## non permuted inputs
depth=200 #seabed depth
fairR=42.5 #fairlead radius
fair_depth=(15) #depth of fairlead

## Permuted inputs
D=np.arange(0.06,0.18,0.06) # lines diameters
angle1=np.arange(0,120,10) # first line heading
angle2=np.arange(120,240,10) # second line heading
angle3=np.arange(240,360,10) # third line heading
ancr=np.arange(240*3,240*5+240,240) # anchor radius
L_ratio=(0.5,0.7,0.9) # L= Lmin+L_ratio(Lmax-Lmin)

## Wind thrust per wind speed
Thrust_wind=pd.read_csv(r'D:\Moordyn\moorpy\MoorPy\Activefloat_AeroThrust.csv',header=0,names=['WindSpeed','Thrust'])
Thrust_wind.sort_values(by=['Thrust'],ignore_index=True,inplace=True,ascending=False)
path0= r"D:\Moordyn\moorpy\MoorPy\Outputs\Github_trials"
Const_dict={
    "Yaw":{'X0': [-10,10],'X': [-10,10],'X_diff': [-10,10]},
            "Roll":{'X0': [-99999999999,99999999999],'X': [-9999999999,9999999999],'X_diff': [-2,2]}
            # "Pitch":{'X0': [-99999999999,99999999999],'X': [-9999999999,9999999999],'X_diff': [-10,10]}
            }        
resultfile='\FinalPermutationResult'
forceresultfile='\ForceResult' #### calculating and printing the forces should be optional
permutationmatrixfile='\PermutationMatrix'

# WindSpeed=10
########## Inputs End ##########


if not os.path.exists(path0):
    os.makedirs(path0)

if not os.path.isfile(path0+permutationmatrixfile+'.csv'):
    PermMatrix=PermutationMatrix(depth, fairR, D, angle1, angle2, angle3, ancr,L_ratio )
    PermMatrix.to_csv(path0+ permutationmatrixfile +'.csv')
else:
    PermMatrix=pd.read_csv(path0+ permutationmatrixfile +'.csv',header=0,names=['PermNum','ang1', 'D1', 'ancR1','L1','ang2', 'D2', 'ancR2','L2','ang3', 'D3', 'ancR3','L3'])


Trymatrix=PermMatrix

# for ws in range(len(Thrust_wind)): ########## This for loop I want to cancel. I want to start only with the wind speed with highest thrust
path= path0+'\WindSpeed_'+str(Thrust_wind.WindSpeed[0])
print('Wind Speed ')
print(Thrust_wind.WindSpeed[0])
if not os.path.exists(path):
    os.makedirs(path)
Faero=(Thrust_wind.Thrust[0],0,0,0,0,0)
aero_angle=np.arange(0,360,5)
F=FAero_windD(Faero,aero_angle)

# loop over wind directions
for j in range(len(aero_angle)): 
    print('Wind Direction ')
    print(aero_angle[j])
    Forces=F[aero_angle[j]]
    
    X0=np.full((6,len(Trymatrix.D1)),0,dtype=float)
    X=np.full((6,len(Trymatrix.D1)),0,dtype=float)
    X_diff=np.full((6,len(Trymatrix.D1)),0,dtype=float)
    dropped=np.full((len(Trymatrix.D1),3),np.nan)
    Fairlead0F=np.full((len(Trymatrix.D1),9),0,dtype=float)
    FairleadF=np.full((len(Trymatrix.D1),9),0,dtype=float)

    # loop over mooring system designs
    for i in range(len(Trymatrix.D1)):
        # pdb.set_trace()
        type_string=("Chain1","Chain2","Chain3") #naming of the mooring type (can be anything)
        LineD=(1.8*Trymatrix.D1[i],1.8*Trymatrix.D2[i],1.8*Trymatrix.D3[i]) #mooring lines equivalent diameter
        dryMass_L=(19.9*1000*Trymatrix.D1[i]**2,19.9*1000*Trymatrix.D2[i]**2,19.9*1000*Trymatrix.D3[i]**2) #dry mass of lines
        EA=(0.854*10**11*Trymatrix.D1[i]**2,0.854*10**11*Trymatrix.D2[i]**2,0.854*10**11*Trymatrix.D3[i]**2) #axial stiffness
        angle=(Trymatrix.ang1[i],Trymatrix.ang2[i],Trymatrix.ang3[i]) #anchor and fairlead angle (maybe change this in the future to add more DOF) #permuted
        anchorR=(Trymatrix.ancR1[i],Trymatrix.ancR2[i],Trymatrix.ancR3[i]) #anchor radius #permuted
        LineLength=(Trymatrix.L1[i],Trymatrix.L2[i],Trymatrix.L3[i]) #lines length #function of depth and anchor radius set minimum and maximum limits

    # Create the mooring line system
        S=makeSampleSystem(depth, type_string, LineD, dryMass_L, EA, angle, anchorR, fair_depth, fairR, LineLength, mass, volume,rCG, AWP, rM)


    ########### Get equilibrium in absence of wind
    ########### This part of getting equilibrium in absence of wind should be done outside the loop over wind directions
    ########### This slows down the code at the moment
    ########### We should switch loops between the wind direction and the design. We should do each design for all wind direction not vice versa
    ###########
    ########### There is repetition for the process for static equilibrium and for equilibrium after applying the aerodynamic force.
        S.initialize() 
        if not S.solveEquilibrium3(no_fail=True):
            continue
            # pdb.set_trace()
        ############################
# Please get rid of the dropped dataframe
        ############################
        ####### constraint to check that there are no vertical forces on the anchor
        for ll in range(len(S.lineList)):
            if ~(S.lineList[ll].LBot> 0):
                dropped[i,2]=Trymatrix.PermNum[i] # I dont use this dropped dataframe anymore we can remove it.
        if ~np.isnan(dropped[i,2]):
            continue
        
        X0[:,i]=S.getPositions(DOFtype="free")
        # Forces should be only saved if asked by the user
        for ff in range(len(S.lineList)):
            Fairlead0F[i,0+3*ff:3+3*ff]=S.lineList[ff].fB
        ########### end Get equilibrium in absence of wind
        
    # Apply the Aerodynamic forces at 8m/s and recalculate equilibrium and repeating again
        S.bodyList[0].f6Ext=Forces
        if not S.solveEquilibrium3(no_fail=True):
            continue
            # pdb.set_trace()
        ############################

        ############################
        ####### constraint to check that there are no vertical forces on the anchor
        for ll in range(len(S.lineList)):
            if ~(S.lineList[ll].LBot> 0):
                dropped[i,2]=Trymatrix.PermNum[i] # I dont use this dropped dataframe anymore we can remove it.
        if ~np.isnan(dropped[i,2]):
            continue
        X[:,i]=S.getPositions(DOFtype="free")
        X_diff[:,i]=X[:,i]-X0[:,i]
        
        # Forces should be only saved if asked by the user
        for ff in range(len(S.lineList)):
            FairleadF[i,0+3*ff:3+3*ff]=S.lineList[ff].fB
        # pdb.set_trace()
    if X0.shape[1] != X.shape[1] != X_diff.shape[1] != i:
        raise ValueError('The position vectors have the wrong dimensions.')
        
    #  Saving the positions and the forces in dictionaries
    # I believe the way I deal with data here is terrible, there should be a better way
    # Should create a function for this. Also again the forces should be optional
    Position={"Surge":{'X0': X0[0,0:],'X': X[0,0:],'X_diff': X_diff[0,0:]},"Sway":{'X0':X0[1,0:],'X': X[1,0:],'X_diff': X_diff[1,0:]},"Heave":{'X0':X0[2,0:],'X': X[2,0:],'X_diff': X_diff[2,0:]},"Roll":{'X0':X0[3,0:],'X': X[3,0:],'X_diff': X_diff[3,0:]}, "Pitch":{'X0':X0[4,0:],'X': X[4,0:],'X_diff': X_diff[4,0:]},"Yaw":{'X0':X0[5,0:],'X': X[5,0:],'X_diff': X_diff[5,0:]}}
    FairLForce={"F0L1":{'X': Fairlead0F[0:,0],'Y': Fairlead0F[0:,1],'Z': Fairlead0F[0:,2]},"F0L2":{'X': Fairlead0F[0:,3],'Y': Fairlead0F[0:,4],'Z': Fairlead0F[0:,5]},"F0L3":{'X': Fairlead0F[0:,6],'Y': Fairlead0F[0:,7],'Z': Fairlead0F[0:,8]},"FL1":{'X': FairleadF[0:,0],'Y': FairleadF[0:,1],'Z': FairleadF[0:,2]}, "FL2":{'X': FairleadF[0:,3],'Y': FairleadF[0:,4],'Z': Fairlead0F[0:,5]},"FL3":{'X': FairleadF[0:,6],'Y': FairleadF[0:,7],'Z': FairleadF[0:,8]}}                
    
    Position_Data=pd.DataFrame(data=Position, columns=['Surge','Sway','Heave','Roll','Pitch','Yaw'])
    Fairlead_Data=pd.DataFrame(data=FairLForce, columns=['F0L1','F0L2','F0L3','FL1','FL2','FL3'])
    
    # I save the displacement results before considering the motion constraints
    # We should consider the constraints first to decrease the size of the result files
    unnested_lst = []
    for col in Position_Data.columns:
        unnested_lst.append(Position_Data[col].apply(pd.Series).stack())
        result = pd.concat(unnested_lst, axis=1, keys=Position_Data.columns)
        result= result.sort_index(axis=0,level=1)
    result.reset_index(inplace=True)
    result.rename(columns={'level_0': 'displacment', 'level_1': 'PermNum','Surge':'Surge','Sway':'Sway','Heave':'Heave','Roll':'Roll','Pitch':'Pitch','Yaw':'Yaw'},inplace=True)
    result.PermNum = pd.DataFrame(np.repeat(Trymatrix.PermNum.values,len(result)/len(Trymatrix),axis=0))
    result.to_csv(path+resultfile+str(aero_angle[j])+'.csv')
    
    # I repeat the same method again for the forces
    # I save the forces results before considering the motion constraints
    # We should consider the constraints first to decrease the size of the result files
    unnested_lst = []
    for col in Fairlead_Data.columns:
        unnested_lst.append(Fairlead_Data[col].apply(pd.Series).stack())
        resultForces = pd.concat(unnested_lst, axis=1, keys=Fairlead_Data.columns)
        resultForces= resultForces.sort_index(axis=0,level=1)
    resultForces.reset_index(inplace=True)
    resultForces.rename(columns={'level_0': 'direction', 'level_1': 'PermNum','F0L1':'F0L1','F0L2':'F0L2','F0L3':'F0L3','FL1':'FL1','FL2':'FL2','FL3':'FL3'},inplace=True)
    resultForces.PermNum = pd.DataFrame(np.repeat(Trymatrix.PermNum.values,len(result)/len(Trymatrix),axis=0))
    resultForces.to_csv(path+forceresultfile+str(aero_angle[j])+'.csv')
    
    # I apply the constraints and drop all designs that dont reach equilibrium or have vertical forces on the anchor
    Trymatrix= MotionConstraints (Const_dict, path, resultfile+str(aero_angle[j])+'.csv',Trymatrix)
    Trymatrix.to_csv(path+ permutationmatrixfile+str(aero_angle[j])+'.csv')
    print('WindDirection done')
    # pdb.set_trace()
print('WindSpeed done')











