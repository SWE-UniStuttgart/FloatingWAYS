# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 14:04:05 2022

@author: mahfouz
"""

import pandas as pd
import numpy as np
import pdb
import floris.tools as wfct
import matplotlib.pyplot as plt
import math

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

def weibull(x, k=2, lam=8.0):
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

def myround(x, base=5):
    return base * round(x/base)

def delta_energy(wd,ws,freq,layout0,newlayoutdict,plot,AEP_initial,delta=1):
    delta_E=pd.DataFrame()
    layout=layout0   
    NewLayoutX_delta= newlayoutdict["X"].copy()
    NewLayoutY_delta= newlayoutdict["Y"].copy()
    NewLayoutX_delta=NewLayoutX_delta.reindex(sorted(NewLayoutX_delta.columns), axis=1)
    NewLayoutY_delta=NewLayoutY_delta.reindex(sorted(NewLayoutY_delta.columns), axis=1)
    
    col=pd.DataFrame(NewLayoutX_delta.columns.astype(str).astype(float),columns=['MoorPy_ang'])
    col['Floris_ang'] = 90-col['MoorPy_ang']+180
    col['Floris_ang'] %= 360
    
    NewLayoutX_delta.set_axis(col['Floris_ang'],axis=1,inplace=True)
    NewLayoutY_delta.set_axis(col['Floris_ang'],axis=1,inplace=True)
    
    NewLayoutX_delta=NewLayoutX_delta.reindex(sorted(NewLayoutX_delta.columns), axis=1)
    NewLayoutY_delta=NewLayoutY_delta.reindex(sorted(NewLayoutY_delta.columns), axis=1)
    NewLayoutX_delta.set_axis(wd,axis=1,inplace=True)
    NewLayoutY_delta.set_axis(wd,axis=1,inplace=True)
    # pdb.set_trace() 
    
    delta_AEPtotal=[]
    total_AEP0=[]
    total_AEP_new=[]
    for ww in range(len(wd)):
        LayoutperDirection=pd.DataFrame(layout).transpose()
        LayoutperDirection.columns=['x0','y0']
        LayoutperDirection['x_new']= LayoutperDirection.x0
        LayoutperDirection['y_new']= LayoutperDirection.y0
        if delta==1:
            LayoutperDirection['x_new']= NewLayoutX_delta.loc[:,wd[ww]]+LayoutperDirection.x0
            LayoutperDirection['y_new']= NewLayoutY_delta.loc[:,wd[ww]]+LayoutperDirection.y0 
        if delta==0:
            LayoutperDirection['x_new']= NewLayoutX_delta.loc[:,wd[ww]]
            LayoutperDirection['y_new']= NewLayoutY_delta.loc[:,wd[ww]]
        fi.reinitialize_flow_field(layout_array=(LayoutperDirection.x0,LayoutperDirection.y0)) 
        AEP0=fi.get_farm_AEP(np.array([wd[ww]]), np.array([ws[ww]]), np.array([freq[ww]]))
        if plot==1:
            hor_plane0 = fi.get_hor_plane(height=fi.floris.farm.turbines[0].hub_height,x_bounds=[-3000,3000],y_bounds=[-3000,3000])
        
        layout_xnew=LayoutperDirection.x_new
        layout_xnew.loc[LayoutperDirection.x_new.isnull()]=LayoutperDirection.x0.loc[LayoutperDirection.x_new.isnull()]
        layout_ynew=LayoutperDirection.y_new
        layout_ynew.loc[LayoutperDirection.y_new.isnull()]=LayoutperDirection.y0.loc[LayoutperDirection.y_new.isnull()]
        fi.reinitialize_flow_field(layout_array=(LayoutperDirection.x_new,LayoutperDirection.y_new)) #layout fillna
        AEP_new=fi.get_farm_AEP(np.array([wd[ww]]), np.array([ws[ww]]), np.array([freq[ww]]))
        if plot==1:
            hor_plane_new = fi.get_hor_plane(height=fi.floris.farm.turbines[0].hub_height,x_bounds=[-3000,3000],y_bounds=[-3000,3000])
        delta_AEP=(AEP_new-AEP0)/(AEP_initial)*100
        delta_AEPtotal.append(delta_AEP)
        total_AEP0.append(AEP0)
        total_AEP_new.append(AEP_new)
        if plot==1:
            mydpi=96
            fig = plt.figure(figsize=(9,6),dpi=mydpi)
            # fig = plt.figure(figsize=(650/mydpi,450/mydpi),dpi=mydpi)
            ax0 = fig.add_subplot(121)
            ax1 = fig.add_subplot(122)
            ax1.set_yticklabels([])
            ax1.set_xlim(-2900,2900)
            ax0.set_xlim(-2900,2900)
            ax1.set_ylim(-2900,2900)
            ax0.set_ylim(-2900,2900)
            ax1.set_xticks(np.arange(-2400,3600,1200))
            ax0.set_xticks(np.arange(-2400,3600,1200))
            ax0.set_yticks(np.arange(-2400,3600,1200))
            ax1.set_yticks(np.arange(-2400,3600,1200))
            ax0.set_xticklabels(np.arange(-2400/240,3600/240,1200/240))
            ax1.set_xticklabels(np.arange(-2400/240,3600/240,1200/240))
            ax0.set_yticklabels(np.arange(-2400/240,3600/240,1200/240))
            wfct.visualization.visualize_cut_plane(hor_plane0, ax=ax0)
            im=wfct.visualization.visualize_cut_plane(hor_plane_new, ax=ax1)
            
            # ax0.plot(LayoutperDirection.x0,LayoutperDirection.y0,'.k')
            ax1.plot(LayoutperDirection.x0,LayoutperDirection.y0,'.k')
            ax0.set_xlabel('x/D', fontdict = {'fontsize' : 12})
            ax0.set_ylabel('y/D', fontdict = {'fontsize' : 12})
            ax0.yaxis.set_label_coords(-0.1, 0.5)
            ax1.set_xlabel('x/D', fontdict = {'fontsize' : 12})
            ax0.set_title('OWFL')
            ax1.set_title('Targeted layout')
            cax = plt.axes([0.92, 0.25, 0.02, 0.5])
            cbar=fig.colorbar(im,cax=cax,shrink=0.9)
            cbar.ax.set_ylabel('Wind speed [m/s]', rotation=270)
            # ax0.legend(
            # ["OWFL"],
            # loc="lower center",
            # bbox_to_anchor=(0.8, 0.85),
            # ncol=1,
            # fontsize=10)
            ax1.legend(
            ["OWFL"],
            loc="lower center",
            bbox_to_anchor=(0.85, 0.88),
            ncol=1,
            fontsize=11)
            plt.subplots_adjust(wspace=0.1)
            plt.savefig(r'D:\Moordyn\moorpy\MoorPy\Windenergy_journal_fig\Targetedlayout_'+'Wind_dir_%s.png'%(wd[ww]),dpi=200)
    delta_E=pd.DataFrame(delta_AEPtotal,columns=['delta_AEP'])
    delta_E['AEP0']=total_AEP0
    delta_E['AEP_new']=total_AEP_new
    delta_E['Wind_dir']=wd
    delta_E['Wind_dir']=90-delta_E['Wind_dir']+180
    delta_E['Wind_dir']%= 360
    delta_E['Wind_dir']=myround(delta_E['Wind_dir'], base=5)
    delta_E.sort_values(by=['delta_AEP'],ignore_index=True,inplace=True)
    return delta_E

def layoutfromMoorSystem(TurbineMoorDesign,Floris_perp,Xdiff,Ydiff): 
    
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

path =r"D:\Moordyn\moorpy\MoorPy\Outputs\trynew\WindSpeed_10"
# Windrose parameters
wd=[0., 22.5, 45., 67.5, 90., 112.5, 135., 157.5, 180., 202.5, 225., 247.5, 270., 292.5, 315., 337.5]
ws=[10, 10,  10,  10,    10,   10,    10,   10,    10,    10,    10,   10,   10,    10,   10,    10]
freq= [.025,  .024,  .029,  .036,.063,  .065,  .100,  .122,.063,  .038,  .039,  .083, .213,  .046,  .032,  .022]

N_circles=2
D=240
spc=5
layout_base=ConcentricCirclesLayout(N_circles,spc,D)   
# wind farm layout same as the one used in Pyoptsparsetry1
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

# wind tuebine parameters and wake model
fi = wfct.floris_interface.FlorisInterface(
    "D:\FLORIS\FLORIS_15MW.json"
)

################ Baseline layout ###############
mydpi=96
fig= plt.figure(figsize=(4.5,4.5),dpi=mydpi)
fontsize = 13
ax = fig.add_subplot()
ax.plot(np.array(layout_base[0])/240, np.array(layout_base[1])/240, "ob")
# ax.plot(layout0[0], layout0[1], ".r")
ax.set_xlabel("x/D", fontsize=fontsize)
ax.set_ylabel("y/D", fontsize=fontsize,loc='center')
# plt.axis("equal")
ax.set_yticks([-10,-5,0,5,10])

ax.xaxis.set_tick_params( labelsize=fontsize)
ax.yaxis.set_tick_params( labelsize=fontsize)
ax.yaxis.set_label_coords(-0.1, 0.5)
plt.grid()
# plt.savefig(r'D:\Moordyn\moorpy\MoorPy\Windenergy_journal_fig\Baseline.png',dpi=200)
plt.show()


################ Optimized layout ###############
mydpi=96
fig= plt.figure(figsize=(5,5),dpi=mydpi)
fontsize = 13
ax = fig.add_subplot()
ax.plot(np.array(layout_base[0])/240, np.array(layout_base[1])/240, "ob")
ax.plot(np.array(layout0[0])/240, np.array(layout0[1])/240, "or")
ax.set_xlabel("x/D", fontsize=fontsize)
ax.set_ylabel("y/D", fontsize=fontsize,loc='center')
# plt.axis("equal")
ax.set_yticks([-10,-5,0,5,10])
ax.xaxis.set_tick_params( labelsize=fontsize)
ax.yaxis.set_tick_params( labelsize=fontsize)
ax.legend(
    ["Baseline", "Optimized"],
    loc="lower center",
    bbox_to_anchor=(0.5, 1.01),
    ncol=2,
    fontsize=fontsize,
        )
ax.yaxis.set_label_coords(-0.1, 0.5)
plt.grid()
plt.savefig(r'D:\Moordyn\moorpy\MoorPy\Windenergy_journal_fig\OWFL.png',dpi=200)
plt.show()


################ Windrose constant wind speed ###############
mydpi=96
legend_kwargs=dict(ncol=1, loc='best', fancybox=True,bbox_to_anchor=(0.4, 0.55, 0.5, 0.5))
fig = plt.figure(figsize=(650/mydpi,450/mydpi),dpi=mydpi)
ax = fig.add_subplot(projection='polar')
# ax.set_yticklabels([])
ax.xaxis.set_tick_params( labelsize=fontsize)
ax.yaxis.set_tick_params( labelsize=fontsize)
ax.set_rticks([0.05,0.1,0.15,0.2,0.25])
ax.set_rlabel_position(65.5)
Rose_onews = wfct.wind_rose.WindRose()
Data_Rose_onews = Rose_onews.make_wind_rose_from_user_dist(np.array(wd), np.array(ws), np.array(freq), wd=np.array(wd),ws=np.arange(3,26,1))
Rose_onews.plot_wind_rose(ws_right_edges=np.array([5, 10, 15, 25]),ax=ax,wd_bins=np.array(wd),legend_kwargs=legend_kwargs)
ax.legend_.remove()
# plt.legend(['10 m/s'],bbox_to_anchor=(0.4, 0.55, 0.5, 0.5))
# plt.savefig(r'D:\Moordyn\moorpy\MoorPy\Windenergy_journal_fig\windroseconstwind.png',dpi=200)
plt.show()


# ############### Windrose all wind speeds ###############
mydpi=96
legend_kwargs=dict(ncol=1, loc='best', fancybox=True,bbox_to_anchor=(0.4, 0.55, 0.5, 0.5))
fig = plt.figure(figsize=(650/mydpi,450/mydpi),dpi=mydpi)
ax = fig.add_subplot(projection='polar')
# ax.set_yticklabels([])
# ax.legend(, ['>10','<10'], **legend_kwargs)
ax.set_rticks([0.05, 0.1, 0.15, 0.2, 0.25])  
# ax.set_yticklabels([0.05, 0.1, 0.15, 0.2, 0.25])
fontsize=12
ax.xaxis.set_tick_params( labelsize=fontsize)
ax.yaxis.set_tick_params( labelsize=fontsize)
ax.set_rlabel_position(65.5)
Windrose_real=make_wind_rose_from_weibull(wd, freq, wd=wd, ws=np.arange(3,26,1))
Rose = wfct.wind_rose.WindRose()
Data_Rose= Rose.make_wind_rose_from_user_dist(np.array(Windrose_real.wd), np.array(Windrose_real.ws),
                                                    np.array( Windrose_real.freq_val), wd=np.array(wd),ws=np.arange(3,26,1))
Rose.plot_wind_rose(ws_right_edges=np.array([10,20]),ax=ax,wd_bins=np.array(wd),legend_kwargs=legend_kwargs)
plt.legend(['> 10 m/s','< 10 m/s'],bbox_to_anchor=(0.4, 0.55, 0.5, 0.5))
plt.savefig(r'D:\Moordyn\moorpy\MoorPy\Windenergy_journal_fig\windrose.png',dpi=200)
plt.show()


############### 
Floris_perp= pd.read_csv(path+ '\OptimumDynamicfarmlayout.csv')
Floris_perp.columns=Floris_perp.columns.astype(str).astype(float)
Floris_perp_xaxis= pd.read_csv(path+ '\OptimumDynamicfarmlayout_xaxis.csv')
Floris_perp_xaxis.columns=Floris_perp.columns.astype(str).astype(float)
Floris_perp_yaxis= pd.read_csv(path+ '\OptimumDynamicfarmlayout_yaxis.csv')
Floris_perp_yaxis.columns=Floris_perp.columns.astype(str).astype(float)
AEP_angle= pd.read_csv(path+ '\AEP_summary_per_winddir.csv') # I think maybe we dont need to read this file at all and get the FLORIS wind dir from the windrose data only
# pdb.set_trace() 
# Input from postprocessing3
MoorPy_disp= pd.read_csv(path+ '\perpdiff.csv').round(decimals=2) # The perpendicular displacement of each mooring design
MooringPermMatrix= pd.read_csv(path+ '\FinalPermMatrix.csv') # the design parameters vof each mooring system
Xdiff_file=path+ '\Xdiff.csv' # The positions achieved by each mooring design in each wind dir (x-axis)
Ydiff_file=path+ '\Ydiff.csv' # The positions achieved by each mooring design in each wind dir (y-axis)
Xdiff=pd.read_csv(Xdiff_file).round(decimals=2)
Ydiff=pd.read_csv(Ydiff_file).round(decimals=2)

fi.reinitialize_flow_field(layout_array=layout0)
AEP_initial=fi.get_farm_AEP(np.array(wd), np.array(ws), np.array(freq))


fi.reinitialize_flow_field(layout_array=layout_base)
AEP_base=fi.get_farm_AEP(np.array(wd), np.array(ws), np.array(freq))

############### Targeted layout ###############
TargetLayoutdict={"X":Floris_perp_xaxis,"Y":Floris_perp_yaxis}
delta_E_target = delta_energy(wd,ws,freq,layout0,TargetLayoutdict,1,AEP_initial,delta=0)    
delta_E_target.sort_values(by=['Wind_dir'],inplace=True, ignore_index=True)
gain0=sum(delta_E_target['delta_AEP'])
print('gain0=')
print(gain0) 


pdb.set_trace() 
############### Final layout ###############
wht=pd.read_csv(path+'/TBMDfinal_new_SNOPT_tol5.csv',index_col=False)
wht=wht.iloc[:,1:4]
Layoutdict=layoutfromMoorSystem(wht,Floris_perp,Xdiff,Ydiff)
Layout_perp=PerpdispfromMoorSystem(wht,Floris_perp,MoorPy_disp)
delta_E_final = delta_energy(wd,ws,freq,layout0,Layoutdict,0,AEP_initial,delta=1) 
delta_E_final.sort_values(by=['Wind_dir'],inplace=True, ignore_index=True)    
gain0=sum(delta_E_final['delta_AEP'])
print('gain0=')
print(gain0)       



delta_E_finalplot=delta_E_final.copy()
delta_E_finalplot['Wind_dir']=90-delta_E_finalplot['Wind_dir']+180+360
delta_E_finalplot['Wind_dir']%=360
delta_E_finalplot.sort_values(by=['Wind_dir'],inplace=True, ignore_index=True)  
delta_E_finalplot['Wind_dir']=wd

delta_E_targetplot=delta_E_target.copy()
delta_E_targetplot['Wind_dir']=90-delta_E_targetplot['Wind_dir']+180+360
delta_E_targetplot['Wind_dir']%=360
delta_E_targetplot.sort_values(by=['Wind_dir'],inplace=True, ignore_index=True)  
delta_E_targetplot['Wind_dir']=wd


############### Final gain ###############
fig = plt.figure(figsize=(650/mydpi,450/mydpi),dpi=mydpi)
ax = fig.add_subplot(projection='polar')
ax.bar(np.deg2rad(delta_E_finalplot['Wind_dir']),delta_E_finalplot['delta_AEP'],width= np.radians(22.5),
       edgecolor="k",color='g')
ax.set_theta_direction(-1)
ax.set_theta_offset(np.pi / 2.0)
ax.set_rticks([0.2,0.4,0.6])
fontsize=12
ax.xaxis.set_tick_params( labelsize=fontsize)
ax.yaxis.set_tick_params( labelsize=fontsize)
ax.set_theta_zero_location("N")
ax.set_xticklabels(["N", "NE", "E", "SE", "S", "SW", "W", "NW"])
ax.set_yticklabels(["0.2", "0.4", "0.6%"])
ax.set_title('Final energy gain',fontsize=14)
ax.set_rlabel_position(65.5)
plt.savefig(r'D:\Moordyn\moorpy\MoorPy\Windenergy_journal_fig\Final_E_gain.png',dpi=200)
plt.show()


############### Target gain ###############
fig = plt.figure(figsize=(650/mydpi,450/mydpi),dpi=mydpi)
ax = fig.add_subplot(projection='polar')
ax.bar(np.deg2rad(delta_E_targetplot['Wind_dir']),delta_E_targetplot['delta_AEP'],width= np.radians(22.5),
       edgecolor="k",color='g')
ax.set_theta_direction(-1)
ax.set_theta_offset(np.pi / 2.0)
ax.set_rticks([0.2,0.4,0.6])
ax.set_theta_zero_location("N")
ax.xaxis.set_tick_params( labelsize=fontsize)
ax.yaxis.set_tick_params( labelsize=fontsize)
ax.set_xticklabels(["N", "NE", "E", "SE", "S", "SW", "W", "NW"])
ax.set_yticklabels(["0.2", "0.4", "0.6%"])
ax.set_title('Targeted energy gain',fontsize=14)
ax.set_rlabel_position(65.5)
plt.savefig(r'D:\Moordyn\moorpy\MoorPy\Windenergy_journal_fig\Targeted_E_gain.png',dpi=200)
plt.show()

############### subplot Target gain + final gain###############

fig = plt.figure(figsize=(9,6),dpi=mydpi)
# fig = plt.figure(figsize=(650/mydpi,450/mydpi),dpi=mydpi)
ax0 = fig.add_subplot(122,projection='polar')
ax1 = fig.add_subplot(121,projection='polar')
ax0.bar(np.deg2rad(delta_E_finalplot['Wind_dir']),delta_E_finalplot['delta_AEP'],width= np.radians(22.5),
       edgecolor="k",color='g')
ax0.set_theta_direction(-1)
ax0.set_theta_offset(np.pi / 2.0)
ax0.set_rticks([0.2,0.4,0.6])
fontsize=12
ax0.xaxis.set_tick_params( labelsize=fontsize)
ax0.yaxis.set_tick_params( labelsize=fontsize)
ax0.set_theta_zero_location("N")
ax0.set_xticklabels(["N", "NE", "E", "SE", "S", "SW", "W", "NW"])
ax0.set_yticklabels(["0.2", "0.4", "0.6%"])
ax0.set_title('Final energy gain',fontsize=14)
ax0.set_rlabel_position(65.5)
ax1.bar(np.deg2rad(delta_E_targetplot['Wind_dir']),delta_E_targetplot['delta_AEP'],width= np.radians(22.5),
       edgecolor="k",color='g')
ax1.set_theta_direction(-1)
ax1.set_theta_offset(np.pi / 2.0)
ax1.set_rticks([0.2,0.4,0.6])
ax1.set_theta_zero_location("N")
ax1.xaxis.set_tick_params( labelsize=fontsize)
ax1.yaxis.set_tick_params( labelsize=fontsize)
ax1.set_xticklabels(["N", "NE", "E", "SE", "S", "SW", "W", "NW"])
ax1.set_yticklabels(["0.2", "0.4", "0.6%"])
ax1.set_title('Targeted energy gain',fontsize=14)
ax1.set_rlabel_position(65.5)
plt.savefig(r'D:\Moordyn\moorpy\MoorPy\Windenergy_journal_fig\Targeted_Final_E_gain.png',dpi=200)
plt.show()



def delta_energyplot(wd,ws,freq,targetlayoutdict,finallayoutdict,plot,AEP_initial,layout0,delta=1):
    delta_E=pd.DataFrame()
    layout=layout0  
    FinalLayoutX_delta= finallayoutdict["X"].copy()
    FinalLayoutY_delta= finallayoutdict["Y"].copy()
    FinalLayoutX_delta=FinalLayoutX_delta.reindex(sorted(FinalLayoutX_delta.columns), axis=1)
    FinalLayoutY_delta=FinalLayoutY_delta.reindex(sorted(FinalLayoutY_delta.columns), axis=1)
    
    col=pd.DataFrame(FinalLayoutX_delta.columns.astype(str).astype(float),columns=['MoorPy_ang'])
    col['Floris_ang'] = 90-col['MoorPy_ang']+180
    col['Floris_ang'] %= 360
    
    FinalLayoutX_delta.set_axis(col['Floris_ang'],axis=1,inplace=True)
    FinalLayoutY_delta.set_axis(col['Floris_ang'],axis=1,inplace=True)
    
    FinalLayoutX_delta=FinalLayoutX_delta.reindex(sorted(FinalLayoutX_delta.columns), axis=1)
    FinalLayoutY_delta=FinalLayoutY_delta.reindex(sorted(FinalLayoutY_delta.columns), axis=1)
    FinalLayoutX_delta.set_axis(wd,axis=1,inplace=True)
    FinalLayoutY_delta.set_axis(wd,axis=1,inplace=True)
    # pdb.set_trace() 
    
    
    TargetLayoutX_delta= targetlayoutdict["X"].copy()
    TargetLayoutY_delta= targetlayoutdict["Y"].copy()
    TargetLayoutX_delta=TargetLayoutX_delta.reindex(sorted(TargetLayoutX_delta.columns), axis=1)
    TargetLayoutY_delta=TargetLayoutY_delta.reindex(sorted(TargetLayoutY_delta.columns), axis=1)
    
    col=pd.DataFrame(TargetLayoutX_delta.columns.astype(str).astype(float),columns=['MoorPy_ang'])
    col['Floris_ang'] = 90-col['MoorPy_ang']+180
    col['Floris_ang'] %= 360
    
    TargetLayoutX_delta.set_axis(col['Floris_ang'],axis=1,inplace=True)
    TargetLayoutY_delta.set_axis(col['Floris_ang'],axis=1,inplace=True)
    
    TargetLayoutX_delta=TargetLayoutX_delta.reindex(sorted(TargetLayoutX_delta.columns), axis=1)
    TargetLayoutY_delta=TargetLayoutY_delta.reindex(sorted(TargetLayoutY_delta.columns), axis=1)
    TargetLayoutX_delta.set_axis(wd,axis=1,inplace=True)
    TargetLayoutY_delta.set_axis(wd,axis=1,inplace=True)
    
    
    delta_AEPtotal=[]
    total_AEP0=[]
    total_AEP_new=[]
    for ww in range(len(wd)):
        LayoutperDirection=pd.DataFrame([TargetLayoutX_delta.loc[:,wd[ww]], TargetLayoutY_delta.loc[:,wd[ww]]]).transpose()
        # pdb.set_trace() 
        LayoutperDirection.columns=['x0','y0']
        LayoutperDirection['x_new']= LayoutperDirection.x0
        LayoutperDirection['y_new']= LayoutperDirection.y0
        
        LayoutperDirection['x_new']= FinalLayoutX_delta.loc[:,wd[ww]]+layout[0]
        LayoutperDirection['y_new']= FinalLayoutY_delta.loc[:,wd[ww]]+layout[1]
        
        fi.reinitialize_flow_field(layout_array=(LayoutperDirection.x0,LayoutperDirection.y0)) 
        AEP0=fi.get_farm_AEP(np.array([wd[ww]]), np.array([ws[ww]]), np.array([freq[ww]]))
        if plot==1:
            hor_plane0 = fi.get_hor_plane(height=fi.floris.farm.turbines[0].hub_height,x_bounds=[-3000,3000],y_bounds=[-3000,3000])
        
        layout_xnew=LayoutperDirection.x_new
        layout_xnew.loc[LayoutperDirection.x_new.isnull()]=LayoutperDirection.x0.loc[LayoutperDirection.x_new.isnull()]
        layout_ynew=LayoutperDirection.y_new
        layout_ynew.loc[LayoutperDirection.y_new.isnull()]=LayoutperDirection.y0.loc[LayoutperDirection.y_new.isnull()]
        fi.reinitialize_flow_field(layout_array=(LayoutperDirection.x_new,LayoutperDirection.y_new)) #layout fillna
        AEP_new=fi.get_farm_AEP(np.array([wd[ww]]), np.array([ws[ww]]), np.array([freq[ww]]))
        if plot==1:
            hor_plane_new = fi.get_hor_plane(height=fi.floris.farm.turbines[0].hub_height,x_bounds=[-3000,3000],y_bounds=[-3000,3000])
        delta_AEP=(AEP_new-AEP0)/(AEP_initial)*100
        delta_AEPtotal.append(delta_AEP)
        total_AEP0.append(AEP0)
        total_AEP_new.append(AEP_new)
        if plot==1:
            mydpi=96
            fig = plt.figure(figsize=(9,6),dpi=mydpi)
            # fig = plt.figure(figsize=(650/mydpi,450/mydpi),dpi=mydpi)
            ax0 = fig.add_subplot(121)
            ax1 = fig.add_subplot(122)
            ax1.set_yticklabels([])
            ax1.set_xlim(-2900,2900)
            ax0.set_xlim(-2900,2900)
            ax1.set_ylim(-2900,2900)
            ax0.set_ylim(-2900,2900)
            ax1.set_xticks(np.arange(-2400,3600,1200))
            ax0.set_xticks(np.arange(-2400,3600,1200))
            ax0.set_yticks(np.arange(-2400,3600,1200))
            ax1.set_yticks(np.arange(-2400,3600,1200))
            ax0.set_xticklabels(np.arange(-2400/240,3600/240,1200/240))
            ax1.set_xticklabels(np.arange(-2400/240,3600/240,1200/240))
            ax0.set_yticklabels(np.arange(-2400/240,3600/240,1200/240))
            wfct.visualization.visualize_cut_plane(hor_plane0, ax=ax0)
            im=wfct.visualization.visualize_cut_plane(hor_plane_new, ax=ax1)
            
            # ax0.plot(LayoutperDirection.x0,LayoutperDirection.y0,'.k')
            ax1.plot(LayoutperDirection.x0,LayoutperDirection.y0,'.k')
            ax0.set_xlabel('x/D', fontdict = {'fontsize' : 12})
            ax0.set_ylabel('y/D', fontdict = {'fontsize' : 12})
            ax0.yaxis.set_label_coords(-0.1, 0.5)
            ax1.set_xlabel('x/D', fontdict = {'fontsize' : 12})
            ax0.set_title('Targeted layout')
            ax1.set_title('Final layout')
            # ax0.legend(
            # ["OWFL"],
            # loc="lower center",
            # bbox_to_anchor=(0.8, 0.85),
            # ncol=1,
            # fontsize=10)
            cax = plt.axes([0.92, 0.25, 0.01, 0.5])
            cbar=fig.colorbar(im,cax=cax,shrink=0.9)
            cbar.ax.set_ylabel('Wind speed [m/s]', rotation=270)
            ax1.legend(
            ["Targeted"],
            loc="lower center",
            bbox_to_anchor=(0.85, 0.88),
            ncol=1,
            fontsize=11)
            plt.subplots_adjust(wspace=0.1)
            plt.savefig(r'D:\Moordyn\moorpy\MoorPy\Windenergy_journal_fig\Target_vs_Finallayout_'+'Wind_dir_%s.png'%(wd[ww]),dpi=200)
    delta_E=pd.DataFrame(delta_AEPtotal,columns=['delta_AEP'])
    delta_E['AEP0']=total_AEP0
    delta_E['AEP_new']=total_AEP_new
    delta_E['Wind_dir']=wd
    delta_E['Wind_dir']=90-delta_E['Wind_dir']+180
    delta_E['Wind_dir']%= 360
    delta_E['Wind_dir']=myround(delta_E['Wind_dir'], base=5)
    delta_E.sort_values(by=['delta_AEP'],ignore_index=True,inplace=True)
    return delta_E



delta_energyplot(wd,ws,freq,TargetLayoutdict,Layoutdict,1,AEP_initial,layout0,delta=1)


mydpi=96
fig = plt.figure(figsize=(4.5,5),dpi=mydpi)
ax = fig.add_subplot(projection='polar')
ax.plot(np.append(np.arange(0,2*np.pi,np.pi/36),0)
        ,np.append(np.abs(MoorPy_disp.loc[3000,:]),np.abs(MoorPy_disp.loc[3000,:])[0]),'g')
rnd=np.random.randint(110,size=(18))
ax.plot(np.append(np.arange(0,2*np.pi,np.pi/9),0)
        ,np.append(rnd,rnd[0]),'r')
# ax.set_theta_direction(-1)
# ax.set_theta_offset(np.pi / 2.0)
fontsize=12
ax.set_rticks([40,80,120])
ax.set_yticklabels([])
ax.set_theta_zero_location("N")
ax.xaxis.set_tick_params( labelsize=fontsize)
ax.yaxis.set_tick_params( labelsize=fontsize)
ax.set_xticklabels(["N", "NE", "E", "SE", "S", "SW", "W", "NW"])
ax.legend(
    ["$\Delta X_p$", "$\Delta X_{tr}$"],
    loc="lower center",
    bbox_to_anchor=(0.5, 1.08),
    ncol=2,
    fontsize=fontsize,
        )
# ax.set_title('Targeted energy gain',fontsize=14)
# ax.set_rlabel_position(65.5)
plt.savefig(r'D:\Moordyn\moorpy\MoorPy\Windenergy_journal_fig\Xp_Xtr.png',dpi=200)
plt.show()



mydpi=96
fig= plt.figure(figsize=(10,5),dpi=mydpi)
fontsize = 13
ax0 = fig.add_subplot(121)
ax0.plot(np.array(layout_base[0])/240, np.array(layout_base[1])/240, "ob")
# ax.plot(layout0[0], layout0[1], ".r")
ax0.set_xlabel("x/D", fontsize=fontsize)
ax0.set_ylabel("y/D", fontsize=fontsize,loc='center')
# plt.axis("equal")
ax0.set_yticks([-10,-5,0,5,10])

ax0.xaxis.set_tick_params( labelsize=fontsize)
ax0.yaxis.set_tick_params( labelsize=fontsize)
ax0.yaxis.set_label_coords(-0.1, 0.5)
plt.grid()

legend_kwargs=dict(ncol=1, loc='best', fancybox=True,bbox_to_anchor=(0.4, 0.55, 0.5, 0.5))
# fig = plt.figure(figsize=(650/mydpi,450/mydpi),dpi=mydpi)
ax1 = fig.add_subplot(122,projection='polar')
# ax.set_yticklabels([])
ax1.xaxis.set_tick_params( labelsize=fontsize)
ax1.yaxis.set_tick_params( labelsize=fontsize)
ax1.set_rticks([0.05,0.1,0.15,0.2])
ax1.set_rlabel_position(65.5)
Rose_onews = wfct.wind_rose.WindRose()
Data_Rose_onews = Rose_onews.make_wind_rose_from_user_dist(np.array(wd), np.array(ws), np.array(freq), wd=np.array(wd),ws=np.arange(3,26,1))
Rose_onews.plot_wind_rose(ws_right_edges=np.array([5, 10, 15, 25]),ax=ax1,wd_bins=np.array(wd),legend_kwargs=legend_kwargs)
ax1.legend_.remove()
plt.savefig(r'D:\Moordyn\moorpy\MoorPy\Windenergy_journal_fig\Baseline_windrose_together.png',dpi=200)
plt.show()



