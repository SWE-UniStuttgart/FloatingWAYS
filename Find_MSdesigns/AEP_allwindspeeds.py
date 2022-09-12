# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 18:16:22 2022

@author: mahfouz
"""
import pandas as pd
import numpy as np
import pdb
import floris.tools as wfct
import matplotlib.pyplot as plt

def myround(x, base=5):
    return base * round(x/base)
def delta_energy(wd,ws,freq,layout0,newlayoutdict,plot,AEP_initial):
    delta_E=pd.DataFrame()
    layout=layout0   
    NewLayoutX_delta= newlayoutdict["X"]
    NewLayoutY_delta= newlayoutdict["Y"]
    
    
    # fi.reinitialize_flow_field(layout_array=layout)
    # AEP_initial=fi.get_farm_AEP(np.array(wd), np.array(ws), np.array(freq))
    
    col=pd.DataFrame(NewLayoutX_delta.columns.astype(str).astype(float),columns=['MoorPy_ang'])
    col['Floris_ang'] = 90-col['MoorPy_ang']+180
    col['Floris_ang'] %= 360
    
    NewLayoutX_delta.set_axis(col['Floris_ang'],axis=1,inplace=True)
    NewLayoutY_delta.set_axis(col['Floris_ang'],axis=1,inplace=True)
    
    
    NewLayoutX_delta=NewLayoutX_delta.reindex(sorted(NewLayoutX_delta.columns), axis=1)
    NewLayoutY_delta=NewLayoutY_delta.reindex(sorted(NewLayoutY_delta.columns), axis=1)
    NewLayoutX_delta.set_axis(wd,axis=1,inplace=True)
    NewLayoutY_delta.set_axis(wd,axis=1,inplace=True)
    
    
    # direction_dict=[]
    # AEP0_total=[]
    # AEP_new_total=[]
    delta_AEPtotal=[]
    # delta_Pwrtotal=pd.DataFrame()
    for ww in range(len(wd)):
        LayoutperDirection=pd.DataFrame(layout).transpose()
        LayoutperDirection.columns=['x0','y0']
        LayoutperDirection['x_new']= LayoutperDirection.x0
        LayoutperDirection['y_new']= LayoutperDirection.y0
        LayoutperDirection['x_new']= NewLayoutX_delta.loc[:,wd[ww]]+LayoutperDirection.x0
        LayoutperDirection['y_new']= NewLayoutY_delta.loc[:,wd[ww]]+LayoutperDirection.y0    
        fi.reinitialize_flow_field(layout_array=(LayoutperDirection.x0,LayoutperDirection.y0)) 
        AEP0=fi.get_farm_AEP(np.array([wd[ww]]), np.array([ws[ww]]), np.array([freq[ww]]))
        if plot==1:
            hor_plane0 = fi.get_hor_plane(height=fi.floris.farm.turbines[0].hub_height)
        # Pwr0=fi.get_turbine_power()
        
        layout_xnew=LayoutperDirection.x_new
        layout_xnew.loc[LayoutperDirection.x_new.isnull()]=LayoutperDirection.x0.loc[LayoutperDirection.x_new.isnull()]
        layout_ynew=LayoutperDirection.y_new
        layout_ynew.loc[LayoutperDirection.y_new.isnull()]=LayoutperDirection.y0.loc[LayoutperDirection.y_new.isnull()]
        fi.reinitialize_flow_field(layout_array=(LayoutperDirection.x_new,LayoutperDirection.y_new)) #layout fillna
        AEP_new=fi.get_farm_AEP(np.array([wd[ww]]), np.array([ws[ww]]), np.array([freq[ww]]))
        if plot==1:
            hor_plane_new = fi.get_hor_plane(height=fi.floris.farm.turbines[0].hub_height)
        # Pwr_new=fi.get_turbine_power()
        delta_AEP=(AEP_new-AEP0)/(AEP_initial)*100
        # delta_Pwr=(np.array(Pwr_new)-np.array(Pwr0))*freq[ww]*8760/(AEP_initial)*100
        # dictionary={
        #     "wd":wd[ww], "Layout":LayoutperDirection, "AEP0":AEP0, "AEP_new":AEP_new, "delta_AEP":delta_AEP
        #             }
        # AEP0_to tal.append(AEP0)
        # AEP_new_total.append(AEP_new)
        # direction_dict.append(dictionary)
        delta_AEPtotal.append(delta_AEP)
        # delta_Pwrtotal.loc[:,wd[ww]]=delta_Pwr
        if plot==1:
            mydpi=96
            fig = plt.figure(figsize=(650/mydpi,450/mydpi),dpi=mydpi)
            ax0 = fig.add_subplot(121)
            ax1 = fig.add_subplot(122)
            wfct.visualization.visualize_cut_plane(hor_plane0, ax=ax0)
            wfct.visualization.visualize_cut_plane(hor_plane_new, ax=ax1)
            ax0.plot(LayoutperDirection.x0,LayoutperDirection.y0,'.k')
            ax1.plot(LayoutperDirection.x0,LayoutperDirection.y0,'.k')
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

def PerpdispfromMoorSystem(TurbineMoorDesign,Floris_perp,MoorPy_perp):     
    TurbineLayout_perp=pd.DataFrame()
    NewLayout_perp=pd.DataFrame([],columns=Floris_perp.columns)
    for i in range(len(Floris_perp)):
        TurbineLayout_perp.loc[0,MoorPy_perp.columns.values]=MoorPy_perp.loc[TurbineMoorDesign.Design[i],:].transpose().values
        col=pd.DataFrame(TurbineLayout_perp.columns.astype(str).astype(float)+TurbineMoorDesign['deltaAngle'][i],columns=['Floris_ang'])
        col['Floris_ang'] %= 360
        col['MoorPy_ang'] = TurbineLayout_perp.columns
        TurbineLayout_perp.rename(columns=col.set_index('MoorPy_ang')['Floris_ang'],inplace=True)
        TurbineLayout_perp=TurbineLayout_perp.reindex(columns=Floris_perp.columns)
        NewLayout_perp.loc[i,:]=TurbineLayout_perp.values
        TurbineLayout_perp=pd.DataFrame()
    return NewLayout_perp

def weibull(x, k=2.5, lam=8.0):
        """
        This method returns a Weibull distribution corresponding to the input
        data array (typically wind speed) using the specified Weibull
        parameters.
        Args:
            x (np.array): List of input data (typically binned wind speed
                observations).
            k (float, optional): Weibull shape parameter. Defaults to 2.5.
            lam (float, optional): Weibull scale parameter. Defaults to 8.0.
        Returns:
            np.array: Weibull distribution probabilities corresponding to
            values in the input array.
        """
        return (k / lam) * (x / lam) ** (k - 1) * np.exp(-((x / lam) ** k))
    
    
def make_wind_rose_from_weibull(winddir,freqperdir, wd=np.arange(0, 360, 5.0), ws=np.arange(0, 26, 1.0),weibull_k=2,weibull_lam=8):
        """
        Populate WindRose object with an example wind rose with wind speed
        frequencies given by a Weibull distribution. The wind direction
        frequencies are initialized according to an example distribution.
        Args:
            wd (np.array, optional): Wind direciton bin centers (deg). Defaults
            to np.arange(0, 360, 5.).
            ws (np.array, optional): Wind speed bin centers (m/s). Defaults to
                np.arange(0, 26, 1.).
        Returns:
            pandas.DataFrame: Wind rose DataFrame containing at least the
            following columns:
                - **wd** (*float*) - Wind direction bin center values (deg).
                - **ws** (*float*) - Wind speed bin center values (m/s).
                - **freq_val** (*float*) - The frequency of occurance of the
                  wind conditions in the other columns.
        """
        # Use an assumed wind-direction for dir frequency
        wind_dir =winddir
        freq_dir = freqperdir

        freq_wd = np.interp(wd, wind_dir, freq_dir)
        freq_ws = weibull(ws, k=weibull_k, lam=weibull_lam)

        freq_tot = np.zeros(len(wd) * len(ws))
        wd_tot = np.zeros(len(wd) * len(ws))
        ws_tot = np.zeros(len(wd) * len(ws))

        count = 0
        for i in range(len(wd)):
            for j in range(len(ws)):
                wd_tot[count] = wd[i]
                ws_tot[count] = ws[j]

                freq_tot[count] = freq_wd[i] * freq_ws[j]
                count = count + 1

        # renormalize
        freq_tot = freq_tot / np.sum(freq_tot)

        # Load the wind toolkit data into a dataframe
        df = pd.DataFrame()

        # Start by simply round and wrapping the wind direction and wind speed
        # columns
        df["wd"] = wd_tot
        df["ws"] = ws_tot

        # Now group up
        df["freq_val"] = freq_tot

        # Save the df at this point
        Windrose = df
        return Windrose
    
fi = wfct.floris_interface.FlorisInterface(
    "D:\FLORIS\FLORIS_15MW.json"
)    
cutin_ws=3
cutout_ws=25
wd=[0., 22.5, 45., 67.5, 90., 112.5, 135., 157.5, 180., 202.5, 225., 247.5, 270., 292.5, 315., 337.5]
# ws=[10, 10,  10,  10,    10,   10,    10,   10,    10,    10,    10,   10,   10,    10,   10,    10]
freq= [.025,  .024,  .029,  .036,.063,  .065,  .100,  .122,.063,  .038,  .039,  .083, .213,  .046,  .032,  .022]    
Windrose=make_wind_rose_from_weibull(wd, freq, wd=wd, ws=np.arange(3,26,1))
layout0=(np.array([ 1323.01325897,  1574.74856648,   358.73554806,  -249.8382204 ,
        -820.65523951, -2275.84473277, -1700.46670389,  2279.14612064,
        2013.41941944,   232.56256817,  -522.55777988, -1290.64147123,
       -2004.19384489, -2339.12376202, -1064.76075803, -1505.32914464,
         -17.95959987,   883.86777135,  1953.93960647]).tolist(),
        np.array([ -460.14213891,   647.74029244,   929.28575322, -2387.11064568,
         160.83796541,  -761.9325208 ,  1406.28619724,   186.06833023,
        1306.3327418 ,  1733.70240513,  2252.14241634,  2023.2825038 ,
         423.5012005 ,  -161.59433565, -1052.91668174, -1801.9056981 ,
        -851.31643136, -2146.11404297, -1394.0855575 ]).tolist())



path0 =r"D:\Moordyn\moorpy\MoorPy\Outputs\trynew"
Floris_perp= pd.read_csv(path0+ '\OptimumDynamicfarmlayout.csv')
Floris_perp.columns=Floris_perp.columns.astype(str).astype(float)
# AEP_angle= pd.read_csv(path0+ '\AEP_summary_per_winddir.csv') # I think maybe we dont need to read this file at all and get the FLORIS wind dir from the windrose data only
# MooringPermMatrix= pd.read_csv(path0+ '\FinalPermMatrix.csv') # the design parameters vof each mooring system


fi.reinitialize_flow_field(layout_array=layout0)
AEP_initial=fi.get_farm_AEP(np.array(Windrose.wd.values), np.array(Windrose.ws.values), np.array(Windrose.freq_val.values))
Gain_Summary=pd.DataFrame(columns=['windspeed','gain'])
for i in range(cutin_ws,cutout_ws+1,1):
    TempWindrose=Windrose.loc[Windrose.ws==i,:].copy()
    new_wd=TempWindrose.wd.values
    new_freq=TempWindrose.freq_val.values
    new_ws=TempWindrose.ws.values
    path =r"D:\Moordyn\moorpy\MoorPy\Outputs\trynew\WindSpeed_"+str(i)
    
    # MoorPy_disp= pd.read_csv(path+ '\perpdiff.csv').round(decimals=2) # The perpendicular displacement of each mooring design
    
    Xdiff_file=path+ '\Xdiff.csv' # The positions achieved by each mooring design in each wind dir (x-axis)
    Ydiff_file=path+ '\Ydiff.csv' # The positions achieved by each mooring design in each wind dir (y-axis)
    Xdiff=pd.read_csv(Xdiff_file).round(decimals=2)
    Ydiff=pd.read_csv(Ydiff_file).round(decimals=2)
    
    wht=pd.read_csv(path0+'/TBMDfinal_new_SNOPT_tol5.csv',index_col=False)
    wht=wht.iloc[:,1:4]
    Layoutdict=layoutfromMoorSystem(wht,Floris_perp,Xdiff,Ydiff)
    # Layout_perp=PerpdispfromMoorSystem(wht,Floris_perp,MoorPy_disp)
    delta_E = delta_energy(new_wd,new_ws,new_freq,layout0,Layoutdict,0,AEP_initial)     
    gain0=sum(delta_E['delta_AEP'])
    print('gain=')
    print(gain0)
    # pdb.set_trace()
    Gain_Summary.loc[i-cutin_ws,'windspeed']=new_ws[0]
    Gain_Summary.loc[i-cutin_ws,'gain']=gain0
    
Total_gain=sum(Gain_Summary['gain'])