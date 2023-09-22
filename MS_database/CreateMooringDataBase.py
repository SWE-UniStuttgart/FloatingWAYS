# -*- coding: utf-8 -*-
"""
Created on Tue Sep  5 14:48:04 2023

@author: mahfouz
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 16:36:56 2023

@author: youssef
"""

"""
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

"""

import numpy as np
from itertools import product 
import pandas as pd
import moorpy as mp
import pdb
import os
import time 
import gc


def PermutationMatrix(depth, fairR, fair_depth, D, angle1, angle2, angle3, ancr, L_ratio):
    """
    This function creates the permutation matrix from the inputs.
    The material is assumed to be chains and the values are built in the function and cannot be changed.
    The number of mooring lines is assumed to be three.
    Inputs:
        depth [m]: depth of sea bed.
        fairR [m]: radial distance from the center of the floater to the mooring lines fairlead
        D [m]: numpy array of the permutations values of diameter for all three lines
        angle1 [deg]: numpy array of the permutations values of line 1 heading
        angle2 [deg]: numpy array of the permutations values of line 2 heading
        angle3 [deg]: numpy array of the permutations values of line 3 heading
        ancr [m]: numpy array of the anchor radii permutation values for all three lines
        L_ratio: numpy array of the permutation values of the mooring line length between the maximum and minimum allowable values.
        
    Outputs:
        Permutation matrix of mooring design parameters (panda dataframe)
        
    """
    print('Creating Permutation matrix')
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
    Lmin=np.round_(np.sqrt((df2[0]-fairR)**2+(depth-fair_depth)**2),decimals=1)
    Lmax=df2[0]+depth-fairR-fair_depth
    PermutMatrix['L1']=Lmin+df2[1]*(Lmax-Lmin)
    df2=pd.DataFrame(PermutMatrix.anc2.tolist(), index= PermutMatrix.index)
    PermutMatrix['ancR2']=df2[0]
    Lmin=np.round_(np.sqrt((df2[0]-fairR)**2+(depth-fair_depth)**2),decimals=1)
    Lmax=df2[0]+depth-fairR-fair_depth
    PermutMatrix['L2']=Lmin+df2[1]*(Lmax-Lmin)
    df2=pd.DataFrame(PermutMatrix.anc3.tolist(), index= PermutMatrix.index)
    PermutMatrix['ancR3']=df2[0]
    Lmin=np.round_(np.sqrt((df2[0]-fairR)**2+(depth-fair_depth)**2),decimals=1)
    Lmax=df2[0]+depth-fairR-fair_depth
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

def makeSampleSystem(depth, type_string, LineD, dryMass_L, EA, angle, anchorR, fair_depth, fairR, LineLength,mass,volume,rCG,AWP,rM):
    '''
    This function creates the mooring system from the mooring design inputs
    The material is assumed to be chains and the values are built in the function and cannot be changed.
    The number of mooring lines is assumed to be three.
    Inputs: 
        depth: integer, seabed depth from seawater level positive [m]
        type_string: string, name of each line to define it's characteristic
        LineD: diameter of each line [m]
        dryMass_: dry mass of the line [kg/m]
        EA: axial stiffness [N]
        angle: angle of the fairlead and anchor [deg]
        anchorR: anchor radius [m]
        fair_depth: depth of the fairlead from sea water level [m]
        fairR: fairlead radius [m]
        LineLength: length of each line [m]
        mass: mass of floater
        volume: volume of floater
        rCG: CG if the floater
        AWP: water plane area
        rm: metacenter relative to the floater 
    Outputs:
        Mooring system
    '''
    # Create blank system object
    MooringSystem = mp.System()
    
    # Set the depth of the system to the depth of the input value
    MooringSystem.depth = depth
    
    # For the number of variables in the LineD variable, add each Linetype to the system
    for i in range(len(type_string)):
        MooringSystem.addLineType(type_string[i], LineD[i], dryMass_L[i], EA[i],name=type_string[i])
    
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

def CreateMoorPySystems(PermMatrix,CG):
    '''
    Creates a list of Mooring Systems in moorpy format.
    The material is assumed to be chains and the values are built in the function and cannot be changed.
    The number of mooring lines is assumed to be three.

    Inputs:
        PermMatrix: Permutation matrix dataframe

    Outputs:
        Mooring systms list in moorpy classes

    '''
    ############# Creating system matrix and equilibrium ############# maybe separate in two
    # t4=time.time()
    #we are trying to calculate X0 outside the main loop so th we save the extra 72 times its getting calculated
    
    #initializing a list to store mooring system objects
    print('============================================================')
    print('Creating mooring systems')
    print('============================================================')
    SampleSystemList=[]
    
    # dropped=np.nan
    #for loop for initializing and storing all the mooring system
    for i in range(len(PermMatrix)):
            # pdb.set_trace()
            type_string=("Chain1","Chain2","Chain3") #naming of the mooring type (can be anything)
            LineD=(1.8*PermMatrix.D1[i],1.8*PermMatrix.D2[i],1.8*PermMatrix.D3[i]) #mooring lines equivalent diameter
            dryMass_L=(19.9*1000*PermMatrix.D1[i]**2,19.9*1000*PermMatrix.D2[i]**2,19.9*1000*PermMatrix.D3[i]**2) #dry mass of lines
            EA=(0.854*10**11*PermMatrix.D1[i]**2,0.854*10**11*PermMatrix.D2[i]**2,0.854*10**11*PermMatrix.D3[i]**2) #axial stiffness
            angle=(PermMatrix.ang1[i],PermMatrix.ang2[i],PermMatrix.ang3[i]) #anchor and fairlead angle (maybe change this in the future to add more DOF) #permuted
            anchorR=(PermMatrix.ancR1[i],PermMatrix.ancR2[i],PermMatrix.ancR3[i]) #anchor radius #permuted
            LineLength=(PermMatrix.L1[i],PermMatrix.L2[i],PermMatrix.L3[i]) #lines length #function of depth and anchor radius set minimum and maximum limits
    
        # Create the mooring line system
            S=makeSampleSystem(depth, type_string, LineD, dryMass_L, EA, angle, anchorR, fair_depth, fairR, LineLength, mass, volume,CG, AWP, rM)
            S.initialize()
            SampleSystemList.append(S)
    
    return SampleSystemList

def StaticEquilibrium_with_Force(SampleSystemList,Faero,Maero,PermMatrix,aero_angle):    
    '''
    This function get the static equilibrium position of the FOWT assuming a force is applied to the floater at sea water level.
    The static equilibrium is calculated for different mooring systems designs and the results are stored for each.
    The systems that don't pass equilibrium or have vertical forces on the anchor are not accepted and filtered out.
    No rotation is happening inside this function.
    
    Inputs:
        SampleSystemList: List of mooring systems design to be attached to the floater for static equilibrium calculations.
        F: Force dataframe for the forces at each wind direction
        aero_angle: numpy array with the wind directions 
        PermMatrix: The permutation matrix used to create the SampleSystemlist
        
    Outputs:
        Returns a dictionary with the position of the FOWT at each DOF at each wind direction 
        and the permutation matrix and systems list that passed the equlibrium
        
    '''
    # I should create a fuction to make sure length F and aero_angle is equal
    #  I should create a function to make sure length SampleSystemList is equal to length PermMatrix
    print('============================================================')
    print('Calculating Static quilibrium positions with external forces')
    print('============================================================')
    Surge=pd.DataFrame()
    Sway=pd.DataFrame()
    Heave=pd.DataFrame()
    Roll=pd.DataFrame()
    Pitch=pd.DataFrame()
    Yaw=pd.DataFrame()
    # loop over wind directions
    
    for j in range(len(aero_angle)): 
        print('***************')
        print('Wind Direction ')
        print(aero_angle[j])
        print('***************')
               
        X=np.full((6,len(SampleSystemList)),np.nan,dtype=float)        
        dropped=np.nan 
        CG=np.array([-0.191783*np.cos(np.deg2rad(aero_angle[j])), -0.191783*np.sin(np.deg2rad(aero_angle[j])), -11.1276])
        # loop over mooring systems for calculating X and Xdiff
        for i in range(len(SampleSystemList)):
            # Apply the Aerodynamic forces at 8m/s and recalculate equilibrium and repeating again
            SampleSystemList[i].bodyList[0].rCG=CG
            x=np.full((6),np.nan,dtype=float)
            for numb in range(10):  
                if numb==0:
                    F=FAero_windD(Faero,Maero,aero_angle[j],rollangle=float(0),pitchangle=float(0),yawangle=float(0))
                else:
                    F=FAero_windD(Faero,Maero,aero_angle[j],rollangle=np.rad2deg(x[3]),pitchangle=np.rad2deg(x[4]),yawangle=np.rad2deg(x[5]))
                
                Forces=F[aero_angle[j]]
                SampleSystemList[i].bodyList[0].f6Ext=Forces
                
                equilibrium_flag=SampleSystemList[i].solveEquilibrium(no_fail=True,maxIter=500)
  
                if not equilibrium_flag:
                    PermMatrix.drop(PermMatrix.loc[i,:].name,inplace=True)
                    x=np.full((6),np.nan,dtype=float)     
                    break            
                ############################
                x=SampleSystemList[i].getPositions(DOFtype="free")

                ############################
                ####### constraint to check that there are no vertical forces on the anchor   
            if not equilibrium_flag:
                SampleSystemList[i]=np.nan
                continue
            
            for ll in range(len(SampleSystemList[i].lineList)):
                if ~(SampleSystemList[i].lineList[ll].LBot> 0):
                    dropped=PermMatrix.PermNum[i] # I dont use this dropped dataframe anymore we can remove it.
            if ~np.isnan(dropped):
                dropped=np.nan
                PermMatrix.drop(PermMatrix.loc[i,:].name,inplace=True)
                SampleSystemList[i]=np.nan             
                continue
            X[:,i]=x

            

            dropped=np.nan
    

        Position_Data=pd.DataFrame({"Surge":X[0,:],"Sway": X[1,:],"Heave":X[2,:],"Roll": X[3,:], "Pitch":X[4,:],"Yaw": X[5,:]})
        Surge.loc[:,aero_angle[j]]=Position_Data.Surge       
        Sway.loc[:,aero_angle[j]]=Position_Data.Sway       
        Heave.loc[:,aero_angle[j]]=Position_Data.Heave        
        Roll.loc[:,aero_angle[j]]=Position_Data.Roll*180/np.pi       
        Pitch.loc[:,aero_angle[j]]=Position_Data.Pitch*180/np.pi
        Yaw.loc[:,aero_angle[j]]=Position_Data.Yaw*180/np.pi
        
        SampleSystemList = [x for x in SampleSystemList if pd.isnull(x) == False]
        PermMatrix=PermMatrix.dropna()
        PermMatrix=PermMatrix.reset_index(drop=True)
        
        Surge=Surge.dropna()
        Surge=Surge.reset_index(drop=True)
        Sway=Sway.dropna()
        Sway=Sway.reset_index(drop=True) 
        Heave=Heave.dropna()
        Heave=Heave.reset_index(drop=True)
        Roll=Roll.dropna()
        Roll=Roll.reset_index(drop=True)
        Pitch=Pitch.dropna()
        Pitch=Pitch.reset_index(drop=True)
        Yaw=Yaw.dropna()
        Yaw=Yaw.reset_index(drop=True)

        print('Number of mooring systems = '+str(len(SampleSystemList)))
                
        print('WindDirection done')
    
    print('WindSpeed done')
    dict_X={
        "Surge":Surge.copy(),
        "Sway":Sway.copy(),
        "Heave":Heave.copy(),
        "Roll":Roll.copy(), 
        "Pitch":Pitch.copy(),
        "Yaw":Yaw.copy(),
        'PermMatrix': PermMatrix,
        'SystemList':pd.DataFrame(SampleSystemList)
            }
    return dict_X

def StaticEquilibrium_without_Force(SampleSystemList,PermMatrix):
    '''
    This function get the static equilibrium position of the FOWT assuming in the absence of any force acting on the system.
    The static equilibrium is calculated for different mooring systems designs and the results are stored for each.
    The systems that don't pass equilibrium or have vertical forces on the anchor are not accepted and filtered out.
    
    Inputs:
        SampleSystemList: List of mooring systems design to be attached to the floater for static equilibrium calculations.
        PermMatrix: The permutation matrix used to create the SampleSystemlist
        
    Outputs:
        Returns a dictionary with the position of the FOWT at each DOF
        and the permutation matrix and systems list that passed the equlibrium
        
    '''
    ############# Get equilibrium and position function #############   
    print('===============================================================')
    print('Calculating Static quilibrium positions without external forces')
    print('===============================================================')
    
    X0=np.full((6,len(SampleSystemList)),np.nan,dtype=float) 
    # Fairlead0F=np.full((6,9),0,dtype=float)
    SampleSystemListcopy=[]
    #loop for calculating X0 for the mooring systems which have passed the equilibrium test
    
    dropped=np.nan
    for i in range(len(SampleSystemList)):

        SampleSystemList[i].bodyList[0].rCG=np.array([-0.191783, 0, -11.1276])
        SampleSystemList[i].bodyList[0].f6Ext=np.zeros(6)
        if not SampleSystemList[i].solveEquilibrium(no_fail=True):
            PermMatrix.loc[i,:]=np.nan
            continue

        for ll in range(len(SampleSystemList[i].lineList)):
            if ~(SampleSystemList[i].lineList[ll].LBot> 0):
                dropped=PermMatrix.PermNum[i] # I dont use this dropped dataframe anymore we can remove it.
        if ~np.isnan(dropped):
            PermMatrix.loc[i,:]=np.nan
            continue
        X0[:,i]=SampleSystemList[i].getPositions(DOFtype="free")
        SampleSystemListcopy.append(SampleSystemList[i])
        dropped=np.nan

    X0Dataframe=pd.DataFrame(data=X0.transpose(),columns=['Surge','Sway','Heave','Roll','Pitch','Yaw'])
    X0Dataframe['PermNum']=PermMatrix.PermNum
    PermMatrix.dropna(inplace=True)
    X0Dataframe.dropna(inplace=True)
    PermMatrix=PermMatrix.reset_index(drop=True)
    X0Dataframe=X0Dataframe.reset_index(drop=True)
    X0Dataframe.loc[:,['Pitch','Roll','Yaw']]=X0Dataframe.loc[:,['Pitch','Roll','Yaw']]*180/np.pi
    SampleSystemList=SampleSystemListcopy.copy()
    dict_X0={
        "Surge":X0Dataframe.Surge,
        "Sway":X0Dataframe.Sway,
        "Heave":X0Dataframe.Heave,
        "Roll":X0Dataframe.Roll, 
        "Pitch":X0Dataframe.Pitch,
        "Yaw":X0Dataframe.Yaw,
        'PermMatrix': PermMatrix,
        'SystemList':pd.DataFrame(SampleSystemList)
        }

    return dict_X0
   
def Clean_StaticEquilibrium_dictionaries(dictionary1,dictionary2):
    '''
    This function takes two dictionaries with the same number of keys and makes sure they both have the same designs
    The function needs to have the PermMatrix as one of the values in each dictionary or it will not work
    
    Inputs:
        Two dictionaries with the same keys and one of the keys and values should be the permutation matrix.
        
    Outputs:
        Two dictionaries with the same Mooring system design

    '''
    if len(dictionary1['PermMatrix'])>len(dictionary2['PermMatrix']):
        for n in range(len(dictionary2.keys())):
            key=list(dictionary1.keys())[n]
            if key=='PermMatrix':
                continue
            dictionary1[key]=dictionary1[key][
                dictionary1['PermMatrix']['PermNum'].isin(dictionary2['PermMatrix']['PermNum'])
                ]
            dictionary1[key]=dictionary1[key].reset_index(drop=True)
        
        dictionary1['PermMatrix']=dictionary1['PermMatrix'][
            dictionary1['PermMatrix']['PermNum'].isin(dictionary2['PermMatrix']['PermNum'])
            ]
        # pdb.set_trace()
        dictionary1['PermMatrix']=dictionary1['PermMatrix'].reset_index(drop=True)
    elif len(dictionary1['PermMatrix'])<len(dictionary2['PermMatrix']):
        for n in range(len(dictionary2.keys())):
            key=list(dictionary2.keys())[n]
            if key=='PermMatrix':
                continue
            # pdb.set_trace()
            dictionary2[key]=dictionary2[key][
                dictionary2['PermMatrix']['PermNum'].isin(dictionary1['PermMatrix']['PermNum'])
                ]
            dictionary2[key]=dictionary2[key].reset_index(drop=True)
        
        dictionary2['PermMatrix']=dictionary2['PermMatrix'][
            dictionary2['PermMatrix']['PermNum'].isin(dictionary1['PermMatrix']['PermNum'])
            ]
        dictionary2['PermMatrix']=dictionary2['PermMatrix'].reset_index(drop=True)
    
    return dictionary1, dictionary2

def delta_Displacement(dictionary1, dictionary2, DOF=['Surge', 'Sway','Heave','Roll','Pitch']):
    '''
    This function subtracts the first dictionary from the second dictionary, for the specified keys (DOFs) 
    equation: Second dictionary - First dictionary
    Both dictionaries should share the same keys mentioned in the DOF, and both dictionaries values should have the same length.
    

    Inputs:
        dictionary1: First dictionary
        dictionary2: Second dictionary
    
    Outputs:
        The resulting difference dictionary for the defined DOF

    '''
    # Should make sure the keys exist in both dictionaries
    # Should make sure that the system list and the permmatrix exist in every dictionary
    
    dictionary_diff={
        'PermMatrix': dictionary2['PermMatrix'],
        'SystemList':dictionary2['SystemList']
        }
    for n in range(len(DOF)):
        key=DOF[n]
        dictionary_diff[key]=dictionary2[key].subtract(dictionary1[key],axis=0)
    
    return dictionary_diff

def Rotating_ref_frame(dictionary, DOF=['Surge', 'Sway','Roll','Pitch']):
    '''
    This function rotates the DOFs given to be in the wind direction and perpendicular to the wind direction
    Inputs:
        dictionary: the dictionary to rotate
        DOF: the degress of fredom to be rotated
        
    Outputs:
        The rotated dictionary with the DOF included in the inputs
    
    '''
    # Yaw and Heave shouldnt be in DOF and this should be checked
    dictionary_rotated={
        'PermMatrix': dictionary['PermMatrix'],
        'SystemList':dictionary['SystemList']
        }
    if 'Surge' in DOF:
        df=pd.DataFrame(np.tile(dictionary['Surge'].columns.astype(str).astype(int),(len(dictionary['Surge']),1)), 
                        columns=dictionary['Surge'].columns)
        Surge=dictionary['Surge']*np.cos(np.deg2rad(df))+dictionary['Sway']*np.sin(np.deg2rad(df))
        dictionary_rotated['Surge']=Surge
        
    if 'Sway' in DOF:
        df=pd.DataFrame(np.tile(dictionary['Sway'].columns.astype(str).astype(int),(len(dictionary['Sway']),1)), 
                        columns=dictionary['Sway'].columns)
        Sway=dictionary['Sway']*np.cos(np.deg2rad(df))-dictionary['Surge']*np.sin(np.deg2rad(df))
        dictionary_rotated['Sway']=Sway
        
    if 'Roll' in DOF:
        df=pd.DataFrame(np.tile(dictionary['Roll'].columns.astype(str).astype(int),(len(dictionary['Roll']),1)), 
                        columns=dictionary['Roll'].columns)
        Roll=dictionary['Roll']*np.cos(np.deg2rad(df))+dictionary['Pitch']*np.sin(np.deg2rad(df))
        dictionary_rotated['Roll']=Roll
        
    if 'Pitch' in DOF:
        df=pd.DataFrame(np.tile(Displacements_dict['Pitch'].columns.astype(str).astype(int),(len(Displacements_dict['Pitch']),1)), 
                        columns=Displacements_dict['Pitch'].columns)
        Pitch=dictionary['Pitch']*np.cos(np.deg2rad(df))-dictionary['Roll']*np.sin(np.deg2rad(df))
        dictionary_rotated['Pitch']=Pitch
    
    
    return dictionary_rotated

def MotionConstraints(dictionary,Const_dict):
    '''
    This function applies the constraints on the input dictionary
    
    Inputs:
        dictionary: the input dictionary
        Const_dict: The constraints dictionary
        
    Outputs:
        The dictionary after applying the contsraints
    '''
    boolean=(((dictionary['Surge']<Const_dict['Surge']['max']) & (dictionary['Surge']>Const_dict['Surge']['min'])) &
            ((dictionary['Sway']<Const_dict['Sway']['max']) & (dictionary['Sway']>Const_dict['Sway']['min'])) &
             ((dictionary['Roll']<=Const_dict['Roll']['max']) & (dictionary['Roll']>=Const_dict['Roll']['min'])) &
             ((dictionary['Pitch']<=Const_dict['Pitch']['max']) & (dictionary['Pitch']>=Const_dict['Pitch']['min'])) &
             ((dictionary['Yaw']<=Const_dict['Yaw']['max']) & (dictionary['Yaw']>=Const_dict['Yaw']['min']))
             )
    for n in range(len(dictionary.keys())):
        key=list(dictionary.keys())[n]
        dictionary[key]=dictionary[key].loc[boolean.all(axis=1),:]
        # pdb.set_trace()
    return dictionary

def Savefiles(dictionary,path,filename,keys=['Surge','Sway','Heave','Roll','Pitch','Yaw'],allkeys=0):
    if allkeys ==1:
        keys=keys
        df=pd.DataFrame(dictionary.get(key) for key in keys).transpose()
        df.to_csv(path+ filename +'.csv', index=False)
    else:    
        for i in range(len(keys)):
            df=pd.DataFrame(dictionary.get(keys[i]))
            df.to_csv(path+ filename+ keys[i] +'.csv', index=False)


########## Inputs Start ##########
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

#### Wind inputs
# Thrust_wind=pd.read_csv(r'D:\Moordyn\moorpy\MoorPy\Activefloat_AeroThrust.csv',header=0,names=['WindSpeed','Thrust'])
# Thrust_wind.sort_values(by=['Thrust'],ignore_index=True,inplace=True,ascending=False)
Faero=np.array([2148821.67906326,-6262.2958047442,-15538.6832881499])
Maero=np.array([17214575.6193484,3892317.0989765,2791002.68346595])
aero_angle=np.arange(0,360,5)
WindSpeed=10

## Permuted inputs
D=np.arange(0.12,0.21,0.03) # lines diameters
angle1=np.arange(0,120,10) # first line heading
angle2=np.arange(120,240,10) # second line heading
angle3=np.arange(240,360,10) # third line heading
ancr=np.arange(240*2.5,240*3.5+240,240) # anchor radius
L_ratio=(0.4,0.6,0.8) # L= Lmin+L_ratio(Lmax-Lmin)

split=0
print('================ 12 total splits ============================================')
print(str(split))
print('=============================================================================')
## Wind thrust per wind speed
path0= r"/home/mahfouz/MooringDatabases/MooringDataBase_UpdateLength"
split_path='/split_'+str(split)
path_oufiles = path0+ split_path +'/WindSpeed_'+str(WindSpeed)
permutationmatrixfile='/PermutationMatrix'
X0finalmatrixfile='/X0finalresult'

Const_dict={
    "Surge":{'min': -9999999,'max': 240.1},
    "Sway":{'min': -9999999,'max': 99999999},
    "Roll":{'min': -9999999,'max': 99999999},
    "Pitch":{'min': -9999999,'max': 99999999},
    "Yaw":{'min': -10,'max': 10}
            } 

########## Inputs End ##########


# F=FAero_windD(Faero,aero_angle)

if not os.path.exists(path0):
    os.makedirs(path0)

if not os.path.isfile(path0+permutationmatrixfile+'.csv'):
    PermMatrix_full=PermutationMatrix(depth, fairR, fair_depth, D, angle1, angle2, angle3, ancr,L_ratio )
    PermMatrix_full.to_csv(path0+ permutationmatrixfile +'.csv', index=False)
else:
    PermMatrix_full=pd.read_csv(path0+ permutationmatrixfile +'.csv',header=0,names=['PermNum','ang1', 'D1', 'ancR1','L1','ang2', 'D2', 'ancR2','L2','ang3', 'D3', 'ancR3','L3'])

if not os.path.exists(path_oufiles):
    os.makedirs(path_oufiles)

# PermMatrix_full = pd.read_csv(path_oufiles + '/NOCONSTRAINTS_PermMatrix.csv')     
PermMatrix_full=np.array_split(PermMatrix_full,12)
Trymatrix = PermMatrix_full[split]
Trymatrix = Trymatrix.reset_index(drop=True)
# pdb.set_trace()
del PermMatrix_full
gc.collect()
SystemList = CreateMoorPySystems(Trymatrix,rCG)
print('Number of mooring systems = '+str(len(SystemList)))
Displacements0_dict = StaticEquilibrium_without_Force(SystemList.copy(),Trymatrix.copy())
Savefiles(Displacements0_dict,path_oufiles,'/X0_NOCONSTRAINTS',keys=['Surge','Sway','Heave','Roll','Pitch','Yaw'],allkeys=1)
Savefiles(Displacements0_dict,path_oufiles,'/NOCONSTRAINTS0_',keys=['PermMatrix'],allkeys=0)

Displacements_dict = StaticEquilibrium_with_Force(list(Displacements0_dict['SystemList'][0]).copy(),Faero,Maero,Displacements0_dict['PermMatrix'].copy(),aero_angle)
Savefiles(Displacements_dict,path_oufiles,'/X_NOCONSTRAINTS_',keys=['Surge','Sway','Heave','Roll','Pitch','Yaw'],allkeys=0)
Savefiles(Displacements_dict,path_oufiles,'/NOCONSTRAINTS_',keys=['PermMatrix'],allkeys=0)

Displacements0_dict, Displacements_dict = Clean_StaticEquilibrium_dictionaries(Displacements0_dict,Displacements_dict)
Displacements_diff_dict = delta_Displacement(Displacements0_dict, Displacements_dict, DOF=['Surge', 'Sway','Heave'])

Displacements_diff_rotated_tran_dict = Rotating_ref_frame(Displacements_diff_dict, DOF=['Surge', 'Sway'])
Displacements_rotated_ang_dict = Rotating_ref_frame(Displacements_dict, DOF=['Roll', 'Pitch'])
Displacements_dict_forconstraints = {
    "Surge":Displacements_diff_rotated_tran_dict['Surge'].copy(),
    "Sway":Displacements_diff_rotated_tran_dict['Sway'].copy(),
    "Heave":Displacements_dict['Heave'].copy(),
    "Roll":Displacements_rotated_ang_dict['Roll'].copy(), 
    "Pitch":Displacements_rotated_ang_dict['Pitch'].copy(),
    "Yaw":Displacements_dict['Yaw'].copy(),
    'PermMatrix': Displacements_dict['PermMatrix'].copy(),
    'SystemList':Displacements_dict['SystemList'].copy()
        }

# pdb.set_trace()

Displacements_dict_constrained = MotionConstraints(Displacements_dict_forconstraints , Const_dict)

Displacements_dict_constrained, Displacements_dict = Clean_StaticEquilibrium_dictionaries(Displacements_dict_constrained , Displacements_dict)
Savefiles(Displacements_dict,path_oufiles,'/X_',keys=['Surge','Sway','Heave','Roll','Pitch','Yaw'],allkeys=0)
Savefiles(Displacements_dict,path_oufiles,'/Final_',keys=['PermMatrix'],allkeys=0)

Displacements_dict_constrained, Displacements_diff_dict = Clean_StaticEquilibrium_dictionaries(Displacements_dict_constrained , Displacements_diff_dict)
Savefiles(Displacements_diff_dict,path_oufiles,'/Xdiff_',keys=['Surge','Sway'],allkeys=0)

Displacements_dict_constrained, Displacements0_dict = Clean_StaticEquilibrium_dictionaries(Displacements_dict_constrained , Displacements0_dict)
Savefiles(Displacements0_dict,path_oufiles,'/X0_',keys=['Surge','Sway','Heave','Roll','Pitch','Yaw'],allkeys=1)

Displacements_dict_constrained, Displacements_diff_rotated_tran_dict = Clean_StaticEquilibrium_dictionaries(Displacements_dict_constrained , Displacements_diff_rotated_tran_dict)


Savefiles(Displacements_diff_rotated_tran_dict,path_oufiles,'/perp_diff_',keys=['Sway'],allkeys=0)
# pdb.set_trace()


