#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 18:33:04 2022

@author: youssef
"""
import time
import pandas as pd
# pd.options.mode.chained_assignment = None  # default='warn'
import numpy as np
import pdb
from floris.tools import FlorisInterface
from floris.tools.visualization import visualize_cut_plane
import matplotlib.pyplot as plt

def myround(x, base=5):
    return base * round(x/base)
def delta_energy(wd,ws,freq,layout0,newlayoutdict,plot,AEP_initial,AEP_dir):
    # t0=time.time()
    delta_E=pd.DataFrame()
    layout=layout0   
    NewLayoutX_delta= newlayoutdict["X"].copy()
    NewLayoutY_delta= newlayoutdict["Y"].copy()
    
    
    col=pd.DataFrame(NewLayoutX_delta.columns.astype(str).astype(float),columns=['MoorPy_ang'])
    col['Floris_ang'] = 90-col['MoorPy_ang']+180
    col['Floris_ang'] %= 360
    
    NewLayoutX_delta.set_axis(col['Floris_ang'],axis=1,inplace=True)
    NewLayoutY_delta.set_axis(col['Floris_ang'],axis=1,inplace=True)
    
    
    NewLayoutX_delta=NewLayoutX_delta.reindex(sorted(NewLayoutX_delta.columns), axis=1)
    NewLayoutY_delta=NewLayoutY_delta.reindex(sorted(NewLayoutY_delta.columns), axis=1)
    NewLayoutX_delta.set_axis(wd,axis=1,inplace=True)
    NewLayoutY_delta.set_axis(wd,axis=1,inplace=True)
    

    delta_AEPtotal=[]
    for ww in range(len(wd)):
        LayoutperDirection=pd.DataFrame()
        LayoutperDirection['x_new']= NewLayoutX_delta.loc[:,wd[ww]]+layout[0]
        LayoutperDirection['y_new']= NewLayoutY_delta.loc[:,wd[ww]]+layout[1]  
  
        AEP0 = AEP_dir.loc[0,wd[ww]]

        fi.reinitialize(layout_x=np.array(LayoutperDirection.x_new).tolist(),
                        layout_y=np.array(LayoutperDirection.y_new).tolist(),
                        wind_directions=[wd[ww]],
                        wind_speeds=ws) #layout fillna

        AEP_new=fi.get_farm_AEP([freq[ww]],freq_warn=0)
        delta_AEP=(AEP_new-AEP0)/(AEP_initial)*100
        delta_AEPtotal.append(delta_AEP)

    delta_E=pd.DataFrame(delta_AEPtotal,columns=['delta_AEP'])
    delta_E['Wind_dir']=wd
    delta_E['Wind_dir']=90-delta_E['Wind_dir']+180
    delta_E['Wind_dir']%= 360
    delta_E['Wind_dir']=myround(delta_E['Wind_dir'], base=5)
    delta_E.sort_values(by=['delta_AEP'],ignore_index=True,inplace=True)

    return delta_E

def layoutfromMoorSystem(TurbineMoorDesign,Floris_perp,Xdiff,Ydiff): 
    # Xdiff=pd.read_csv(Xdiff_file).round(decimals=2)
    # Ydiff=pd.read_csv(Ydiff_file).round(decimals=2)
    
    TurbineLayoutX=pd.DataFrame()
    TurbineLayoutY=pd.DataFrame()
    NewLayoutX=pd.DataFrame([],columns=Floris_perp.columns)
    NewLayoutY=pd.DataFrame([],columns=Floris_perp.columns)
    for i in range(len(Floris_perp)):
        TurbineLayoutX.loc[0,Xdiff.columns.values]=Xdiff.loc[TurbineMoorDesign.Design[i],:].transpose().values
        TurbineLayoutY.loc[0,Ydiff.columns.values]=Ydiff.loc[TurbineMoorDesign.Design[i],:].transpose().values
        X0new_Surge=TurbineLayoutX*np.cos(np.deg2rad(-TurbineMoorDesign['deltaAngle'][i]))+TurbineLayoutY*np.sin(np.deg2rad(-TurbineMoorDesign['deltaAngle'][i]))
        X0new_Sway=TurbineLayoutY*np.cos(np.deg2rad(-TurbineMoorDesign['deltaAngle'][i]))-TurbineLayoutX*np.sin(np.deg2rad(-TurbineMoorDesign['deltaAngle'][i]))
        TurbineLayoutX=X0new_Surge
        TurbineLayoutY=X0new_Sway
        # this doesnt change the angle from moorpy to floris it only rotates the mooring angles
        # The naming convention here needs changing
        # This is still the conventional x-y axis not the floris angles convention
        col=pd.DataFrame(TurbineLayoutX.columns.astype(str).astype(float)+TurbineMoorDesign['deltaAngle'][i],columns=['Floris_ang'])
        col['Floris_ang'] %= 360
        col['MoorPy_ang'] = TurbineLayoutX.columns
        TurbineLayoutX.rename(columns=col.set_index('MoorPy_ang')['Floris_ang'],inplace=True)
        TurbineLayoutY.rename(columns=col.set_index('MoorPy_ang')['Floris_ang'],inplace=True)
        TurbineLayoutX=TurbineLayoutX.reindex(columns=Floris_perp.columns)
        TurbineLayoutY=TurbineLayoutY.reindex(columns=Floris_perp.columns)
        NewLayoutX.loc[i,:]=TurbineLayoutX.values
        NewLayoutY.loc[i,:]=TurbineLayoutY.values
        TurbineLayoutX=pd.DataFrame()
        TurbineLayoutY=pd.DataFrame()    
    # NewLayoutX.to_csv(outputlayoutxfile, index=False)
    # NewLayoutY.to_csv(outputlayoutyfile, index=False)
    NewLayoutdict={"X":NewLayoutX,"Y":NewLayoutY}
    return NewLayoutdict
    # print('Layout after considering mooring systems is saved')
    
# path= r"D:\Moordyn\moorpy\MoorPy\Outputs"
# path =r"D:\Moordyn\moorpy\MoorPy\Outputsnewlayout"

########### Inputs start  
mindist = 4
spc = 6 # maximum spacing at the beginning of the simulation for boundary creation
shape ='circle'
no_circles_squares =1 
no_turbines = 7 
file_number = f'{mindist}{mindist}{no_turbines}{no_turbines}'
windrose = 'iea'

Ti = 0.06

############################ 
path0=fr"../Results_{windrose}"
path =fr"../Results_{windrose}/Results_{mindist}_{spc}"
opt_baseline_file ='/SNOPTlayoutnew_'+shape+'_'+str(no_turbines)+'T_iea.csv'
opt_dyn_layout_file='/OptimumDynamicfarmlayout_'+shape+'_'+str(no_turbines)+'T_iea_TI_'+str(Ti)+'.csv'
opt_aep_sum_file ='/AEP_summary_per_winddir_'+shape+'_'+str(no_turbines)+'T_iea_TI_'+str(Ti)+'.csv'

outputfile='/TBMDfinal_Yaw5deg_'+shape+'_'+str(no_turbines)+'T_iea.csv'
outputfile_process='/TBMD_notfinal_Yaw5deg_'+shape+'_'+str(no_turbines)+'T_iea.csv'

# Mooring Database inputs
MoorPy_disp= pd.read_csv(path0+'/MooringDatabase_Yaw5deg/WindSpeed_10'+ '/perp_diff_Sway.csv').round(decimals=2) # The perpendicular displacement of each mooring design
MooringPermMatrix= pd.read_csv(path0+'/MooringDatabase_Yaw5deg/WindSpeed_10'+ '/Final_PermMatrix.csv') # the design parameters vof each mooring system
Xdiff_file=path0+'/MooringDatabase_Yaw5deg/WindSpeed_10'+ '/Xdiff_Surge.csv' # The positions achieved by each mooring design in each wind dir (x-axis)
Ydiff_file=path0+'/MooringDatabase_Yaw5deg/WindSpeed_10'+ '/Xdiff_Sway.csv' # The positions achieved by each mooring design in each wind dir (y-axis)

MooringPermMatrix['ang1']=myround(MooringPermMatrix['ang1'], base=10)
MooringPermMatrix['ang2']=myround(MooringPermMatrix['ang2'], base=10)
MooringPermMatrix['ang3']=myround(MooringPermMatrix['ang3'], base=10)

# pdb.set_trace()

print('================================')
print(f'{windrose} {no_turbines} Turbines Yaw 5 deg')
print('================================')

# Windrose parameters
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

# wind farm layout same as the one used in Pyoptsparsetry1
layoutdf=pd.read_csv(path+ '/BaselineOptimization' +opt_baseline_file)

layout0=(np.array(layoutdf.x).tolist(),
        np.array(layoutdf.y).tolist())

# wind tuebine parameters and wake model
fi = FlorisInterface('../Inputs/gch15MW.yaml')
fi.floris.flow_field.turbulence_intensity=0.06
# # Input from Pyoptsparsetry1
# Floris_X= pd.read_csv(path+ '\Position_X_FLORIS.csv') # the positions we want the turbines to be at every wind direction (x-axis)
# Floris_Y= pd.read_csv(path+ '\Position_Y_FLORIS.csv') # the positions we want the turbines to be at every wind direction (y-axis)
Floris_perp= pd.read_csv(path+ '/SingleWIndDir_Opt'+opt_dyn_layout_file)
AEP_angle= pd.read_csv(path+ '/SingleWIndDir_Opt'+opt_aep_sum_file) # I think maybe we dont need to read this file at all and get the FLORIS wind dir from the windrose data only
# pdb.set_trace() 

# Mooring Database inputs
Xdiff=pd.read_csv(Xdiff_file).round(decimals=2)
Ydiff=pd.read_csv(Ydiff_file).round(decimals=2)

tolerance=5
########### Inputs end 

Floris_perp.columns=Floris_perp.columns.astype(str).astype(float)


########### calculating the initial AEP
# t0=time.time()
fi.reinitialize(layout=layout0,
                wind_directions=wd,
                wind_speeds=ws)
# t1=time.time()
AEP_initial=fi.get_farm_AEP(freq=freq)
# t2=time.time()
AEP_dir=pd.DataFrame()
for ww in range(len(wd)):

    fi.reinitialize(layout_x=layout0[0],
                    layout_y=layout0[1],
                    wind_directions=[wd[ww]],
                    wind_speeds=ws) 
    # t2=time.time()
    
    AEP_dir.loc[0,wd[ww]]=fi.get_farm_AEP([freq[ww]],freq_warn=0)
    
print("=====================================================")
print('AEP_initial='+str(AEP_initial*10**-9))
print("=====================================================")
# ww=0
# t0=time.time()
# fi.reinitialize(layout=layout0,
#                 wind_directions=[wd[ww]],
#                 wind_speeds=ws) 
# t1=time.time()

# AEP0=fi.get_farm_AEP([freq[ww]],freq_warn=0)

# t2=time.time()

TurbinesAngles=pd.DataFrame(np.tile(AEP_angle['MoorPy_angles_rnd'].transpose(),(len(Floris_perp),1))).transpose()
TurbinesAngles.columns=TurbinesAngles.columns.astype(str).astype(float)



summarygain0=pd.DataFrame()
summarygain1=pd.DataFrame()
WinddirTurbinedict={}



###### We should change this loop to be first wind turbines and then wind direction
###### loop over wind direction
for count in range(len(Floris_perp.columns)):

    TMD_designnum=pd.DataFrame()
    TMD_deltaangle=pd.DataFrame()
    Turbinedict={}
    
###### Can be a function on its own
###### This loop looks for every displacement for each turbine at each wind direction and finds the mooring systems that have this value
###### It also finds the angle where this value is achieved
###### so finding the values matching Sungle_disp with MoorPy_disp

###### There might be a better way to do this as there are only 5 possible displacement (-120,-60,0,60,120). There will always be the same 
######      number of mooring systems for each of these displacements
###### However the angle at which they happen will change according to the wind direction
###### At the end we save the design number and the delta angle.
    ###### loop over wind turbines
    for t in range(len(Floris_perp)):
        tol=tolerance
        MooringDesign=pd.DataFrame()
        WindDirections=pd.DataFrame()
    
        Single_disp=Floris_perp.loc[t,TurbinesAngles.iloc[count,0]] # all turbines should have the same order for wind dir
        trryyy=((MoorPy_disp>(Single_disp-tol)) & 
                (MoorPy_disp<(Single_disp+tol)))
        s = trryyy.stack()
        s = s[s].index.to_frame(index=False).set_axis(['idx','cols'], axis=1)
        MooringDesign[TurbinesAngles.iloc[count,0]]=s.idx
        WindDirections[TurbinesAngles.iloc[count,0]]=s.cols
        if len(MooringDesign)==0:
            print('Tolerace is too small')
            # should create an error here and stop the code

        TMD_designnum=pd.concat([TMD_designnum,MooringDesign.loc[:,TurbinesAngles.iloc[count,0]]], ignore_index=True, axis=1)
        saveang=TurbinesAngles.iloc[count,0]-WindDirections[TurbinesAngles.iloc[count,0]].astype(str).astype(float)
        TMD_deltaangle=pd.concat([TMD_deltaangle,saveang], ignore_index=True, axis=1) # The angle the mooring system needs to rotate
###### Can be a function on its own

    
    ######### difference between the distances 
    # We calculate the absolute distance between what the mooring system can achieve and what we demand.
    # The mooring systems are rotated using the delta angle and the difference is calculated
    # Rotation here means we rename the columns with the new angles (old columns + deltaAngles)
    # We output Turbine dataframe with columns ['Design', 'deltaAngle', 'totaldiff']
    # This dataframe will be used throught the rest of the code 
    
    ###### loop over wind turbines
    for t in range(len(Floris_perp)):
        Turbine=pd.DataFrame()
        Turbine['Design']=TMD_designnum[t].copy().dropna()
        Turbine['deltaAngle']=TMD_deltaangle[t].copy().dropna()
        # totaldiff=pd.DataFrame()
        totaldiff=pd.DataFrame()
        # loop over unique values of delta angle
        for j in range(len(Turbine['deltaAngle'].unique())):
            Designs=MoorPy_disp.reindex(Turbine.Design[Turbine.deltaAngle==Turbine['deltaAngle'].unique()[j]]) # finding mooring designs that share the same delta angle
            col=pd.DataFrame(Designs.columns.astype(str).astype(float)+Turbine['deltaAngle'].unique()[j],columns=['Floris_ang']) # new columns names of mooring designs after rotating them with delta angle  
            col['Floris_ang'] %= 360 # getting column names between 0 and 30 degrees
            col['MoorPy_ang'] = Designs.columns # columns names of designs that have same design (the column names are the same for all designs so sure there is a better way to do this)
            Designs.rename(columns=col.set_index('MoorPy_ang')['Floris_ang'],inplace=True) # renaming the columns with the rotated values
            Designs=Designs.reindex(columns=Floris_perp.columns) # reordering columns
            # totaldiff=totaldiff.append(abs(Designs-Floris_perp.loc[t,:])) # difference between what the mooring system can achieve and what we want it to achieve
            totaldiff=pd.concat([totaldiff,abs(Designs-Floris_perp.loc[t,:])]) # difference between what the mooring system can achieve and what we want it to achieve
            # diffX=pd.DataFrame()
            # diffX=diffX.append(abs(Designs-Floris_perp.loc[t,:]))
            diffX=pd.DataFrame()
            diffX=pd.concat([diffX,abs(Designs-Floris_perp.loc[t,:])])
            # pdb.set_trace()
            Turbine.loc[Turbine[(Turbine.Design.isin(totaldiff.index))&(Turbine.deltaAngle==Turbine['deltaAngle'].unique()[j])].index,'totaldiff']=np.sum(diffX,axis=1).values # this line maybe needs to get out of the loop
            
            # pdb.set_trace() 
        # pdb.set_trace() 

        ###### Can be a function on its own
        ###### We pick the design with the minimum distance diff and unique heading angles
        ###### There is an option also to include unique delta angle and headings
        ###### We do this to decrease the number of possible mooring designs
        
        Turbine.sort_values(by=['totaldiff'],ignore_index=True,inplace=True)
        MooringDesigns=MooringPermMatrix.reindex(Turbine.Design)
        MooringDesigns.drop(columns=['PermNum'],inplace=True)
        MooringDesigns['Design']=MooringDesigns.index
        MooringDesigns.index=Turbine.index
        MooringDesigns['deltaAngle']=Turbine['deltaAngle'].copy()
        # pdb.set_trace()
        Turbine=Turbine.iloc[MooringDesigns.drop_duplicates(subset=['ang1','ang2','ang3'],keep='first').index,:]
        # Turbine=Turbine.iloc[MooringDesigns.drop_duplicates(subset=['ang1','ang2','ang3','deltaAngle'],keep='first').index,:] #should be a flag to use this or the line above it
        Turbine=Turbine.reset_index(drop=True)
        ###### Can be a function on its own
        Turbinedict["Turbine_%s" %(t)]=Turbine
    # pdb.set_trace() 

    WinddirTurbinedict["Winddir_%s" %(int(TurbinesAngles.iloc[count,0]))]=Turbinedict
    
    
# pdb.set_trace() 
###### This just reformalates the way WinddirTurbinedict such that it's the turbine and then another dict including the wind direction.
###### In the previous part we should loop first over wind turbines then over wind directions
###### I believe there is a way to directly write Turbinefullsysall without going through writing WinddirTurbinedict
Turbinefullsysall={} 
TBMD_max=pd.DataFrame([],columns=WinddirTurbinedict["Winddir_0"]["Turbine_0"].columns) 
for t in range(len(Floris_perp)):
    Turbinenum_total=[]
    Turbine_i=pd.DataFrame()
    for i in range(len(WinddirTurbinedict)):
        Turbinenum=WinddirTurbinedict["Winddir_%s" %(int(TurbinesAngles.loc[i,0]))]["Turbine_%s" %(t)].copy()
        Turbinenum_total.append(Turbinenum)
        Turbine_i=pd.concat(Turbinenum_total, ignore_index=True)
    # pdb.set_trace()
    Turbine_i['deltaAngle'].loc[Turbine_i['deltaAngle']<0]=Turbine_i['deltaAngle'].loc[Turbine_i['deltaAngle']<0]+360
    # pdb.set_trace()
    Turbine_i=Turbine_i.drop_duplicates(keep='first')
    
    Turbine_i.sort_values(by=['totaldiff'], axis=0,ascending=True,inplace=True,ignore_index=True)    
    Turbinefullsysall["Turbine_%s" %(t)]=Turbine_i.copy()
    TBMD_max.loc[t,:]=Turbinefullsysall["Turbine_%s" %(t)].loc[0,:] # No more need to start TBMD_max with the designs highest totaldiff
    # pdb.set_trace()



# pdb.set_trace()
######## Just to define gain0. No other need for this
# t1=time.time()
Layoutdict=layoutfromMoorSystem(TBMD_max,Floris_perp,Xdiff,Ydiff)
# t2=time.time()
delta_E = delta_energy(wd,ws,freq,layout0,Layoutdict,0,AEP_initial,AEP_dir)  
# t3=time.time()   
gain0=sum(delta_E['delta_AEP'])
print('gain0=')
print(gain0)

# pdb.set_trace()

######## Iterating over mooring designs of each turbine while the others are fixed till we find the maximum gain possible.
######## The while loop to ensure that the gain stopped increasing before we stop the iteration.
###### Can be a function on its own
TBMDfinal=TBMD_max.copy()

###################

### new part only now
# TBMDfinal=pd.read_csv(path+'/FinalMooringWFDesign'+outputfile_process)
# Layoutdict=layoutfromMoorSystem(TBMDfinal,Floris_perp,Xdiff,Ydiff)
# delta_E = delta_energy(wd,ws,freq,layout0,Layoutdict,0,AEP_initial,AEP_dir)  
# gain0=sum(delta_E['delta_AEP'])
# print('gain0=')
# print(gain0)

###################
###################
gainnew=10
gainprev=0

while gainnew != gainprev:       
    gainprev=gain0
    TBMD_temp=TBMDfinal.copy()

    for t in range(len(Turbinefullsysall)):
        
        for dd in range(len(Turbinefullsysall["Turbine_%s" %(t)])):
            TBMD_temp.loc[t,:]=Turbinefullsysall["Turbine_%s" %(t)].loc[dd,:]
            Layoutdict=layoutfromMoorSystem(TBMD_temp,Floris_perp,Xdiff,Ydiff)
            delta_E = delta_energy(wd,ws,freq,layout0,Layoutdict,0,AEP_initial,AEP_dir)
            gain1=sum(delta_E['delta_AEP'])
            
            if gain1>=gain0:
                TBMDfinal=TBMD_temp.copy()
                gain0=gain1
            else:
                TBMD_temp=TBMDfinal.copy()   
        print('gain0=')
        print(gain0)
        TBMDfinal.to_csv(path+'/FinalMooringWFDesign'+outputfile_process, index=False)
    gainnew=gain0 
    print('gainnew=')
    print(gainnew)
    TBMDfinal.to_csv(path+'/FinalMooringWFDesign'+outputfile_process, index=False)

TBMDfinal.to_csv(path+'/FinalMooringWFDesign'+outputfile, index=False)

# pdb.set_trace()

################### function to read output is needed
# wht=pd.read_csv(path+'/TBMDfinal_short_aftercleaning_tol25.csv',index_col=False)
# wht=pd.read_csv(path+'/TBMDfinal_new_SNOPT_tol5.csv',index_col=False)
# wht=wht.iloc[:,1:4]
# Layoutdict=layoutfromMoorSystem(wht,Floris_perp,Xdiff,Ydiff)
# delta_E = delta_energy(wd,ws,freq,layout0,Layoutdict,0,AEP_initial)     
# gain0=sum(delta_E['delta_AEP'])
# print('gain0=')
# print(gain0)       