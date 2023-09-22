# -*- coding: utf-8 -*-
"""
Created on Tue Jul  5 14:04:05 2022

@author: mahfouz
"""
import time
import pandas as pd
import numpy as np
import pdb
from floris.tools import FlorisInterface
from floris.tools import WindRose
from floris.tools.visualization import visualize_cut_plane
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
    
    NewLayoutX_delta=NewLayoutX_delta.set_axis(col['Floris_ang'],axis=1)
    NewLayoutY_delta=NewLayoutY_delta.set_axis(col['Floris_ang'],axis=1)
    
    NewLayoutX_delta=NewLayoutX_delta.reindex(sorted(NewLayoutX_delta.columns), axis=1)
    NewLayoutY_delta=NewLayoutY_delta.reindex(sorted(NewLayoutY_delta.columns), axis=1)
    NewLayoutX_delta=NewLayoutX_delta.set_axis(wd,axis=1)
    NewLayoutY_delta=NewLayoutY_delta.set_axis(wd,axis=1)
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
            # pdb.set_trace() 
            LayoutperDirection['x_new']= NewLayoutX_delta.loc[:,wd[ww]]
            LayoutperDirection['y_new']= NewLayoutY_delta.loc[:,wd[ww]]
        # LayoutperDirection.x0=np.around(LayoutperDirection.x0,0)
        # LayoutperDirection.y0=np.around(LayoutperDirection.y0,0)
        fi.reinitialize(layout_x=np.array(LayoutperDirection.x0).tolist(),
                        layout_y=np.array(LayoutperDirection.y0).tolist(),
                        wind_directions=[wd[ww]],
                        wind_speeds=ws) 
        # pdb.set_trace()
        # t1=time.time()
        AEP0=fi.get_farm_AEP([freq[ww]],freq_warn=0)
        # t2=time.time()
        # print(t2-t1)
        # pdb.set_trace()
        if plot==1:
            hor_plane0 = fi.calculate_horizontal_plane(height=fi.floris.farm.turbine_definitions[0]['hub_height'],x_bounds=[-4000,4000],y_bounds=[-4000,4000])
        # pdb.set_trace() 
        # layout_xnew=LayoutperDirection.x_new
        # layout_xnew.loc[LayoutperDirection.x_new.isnull()]=LayoutperDirection.x0.loc[LayoutperDirection.x_new.isnull()]
        # layout_ynew=LayoutperDirection.y_new
        # layout_ynew.loc[LayoutperDirection.y_new.isnull()]=LayoutperDirection.y0.loc[LayoutperDirection.y_new.isnull()]
        fi.reinitialize(layout_x=np.array(LayoutperDirection.x_new).tolist(),
                        layout_y=np.array(LayoutperDirection.y_new).tolist(),
                        wind_directions=[wd[ww]],
                        wind_speeds=ws) #layout fillna
        AEP_new=fi.get_farm_AEP([freq[ww]],freq_warn=0)
        if plot==1:
            hor_plane_new = fi.calculate_horizontal_plane(height=fi.floris.farm.turbine_definitions[0]['hub_height'],x_bounds=[-4000,4000],y_bounds=[-4000,4000])
            # pdb.set_trace()
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
            ax1.set_xlim(-1920,1920)
            ax0.set_xlim(-1920,1920)
            ax1.set_ylim(-1920,1920)
            ax0.set_ylim(-1920,1920)
            ax1.set_xticks(np.arange(-1440,2800,1440))
            ax0.set_xticks(np.arange(-1440,2800,1440))
            ax0.set_yticks(np.arange(-1440,2800,1440))
            ax1.set_yticks(np.arange(-1440,2800,1440))
            ax0.set_xticklabels(np.arange(-1440/240,2800/240,1440/240))
            ax1.set_xticklabels(np.arange(-1440/240,2800/240,1440/240))
            ax0.set_yticklabels(np.arange(-1440/240,2800/240,1440/240))
            visualize_cut_plane(hor_plane0, ax=ax0)
            # visualize_cut_plane(hor_plane0, ax=ax1)
            im=visualize_cut_plane(hor_plane_new, ax=ax1)
            for i in range(len(LayoutperDirection.x_new)):
                ax1.text(LayoutperDirection.x0[i]+300*np.sin(np.deg2rad(wd[ww])),
                         LayoutperDirection.y0[i]+300*np.cos(np.deg2rad(wd[ww])),
                         f'T {i+1}',
                         fontsize=10)
                ax0.text(LayoutperDirection.x0[i]+300*np.sin(np.deg2rad(wd[ww])),
                         LayoutperDirection.y0[i]+300*np.cos(np.deg2rad(wd[ww])),
                         f'T {i+1}',
                         fontsize=10)
            # ax0.plot(LayoutperDirection.x0,LayoutperDirection.y0,'.k')
            ax1.plot(LayoutperDirection.x0,LayoutperDirection.y0,'ok')
            # ax1.plot(LayoutperDirection.x_new,LayoutperDirection.y_new,'.g')
            ax0.plot(LayoutperDirection.x0,LayoutperDirection.y0,'ok')
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
            bbox_to_anchor=(0.0, 0.92),
            ncol=1,
            fontsize=11)
            plt.subplots_adjust(wspace=0.1)
            # plt.savefig(r'/home/youssef/Documents/Conferences/WESC/Figures/Finallayout_'+'Wind_dir_%s.png'%(wd[ww]),dpi=200)
            plt.savefig(path+'/Figures_WES'+ '/Targetedlayout_'+'Wind_dir_%s.png'%(wd[ww]),dpi=200)
            # plt.savefig(r'D:\Moordyn\moorpy\MoorPy\Windenergy_journal_fig\Targetedlayout_'+'Wind_dir_%s.png'%(wd[ww]),dpi=200)
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
        
        fi.reinitialize(layout_x=np.array(LayoutperDirection.x0).tolist(),
                        layout_y=np.array(LayoutperDirection.y0).tolist(),
                        wind_directions=[wd[ww]],
                        wind_speeds=ws) 
        # pdb.set_trace()

        AEP0=fi.get_farm_AEP([freq[ww]],freq_warn=0)
        if plot==1:
            hor_plane0 = fi.calculate_horizontal_plane(height=fi.floris.farm.turbine_definitions[0]['hub_height'],x_bounds=[-3000,3000],y_bounds=[-3000,3000])
        
        layout_xnew=LayoutperDirection.x_new
        layout_xnew.loc[LayoutperDirection.x_new.isnull()]=LayoutperDirection.x0.loc[LayoutperDirection.x_new.isnull()]
        layout_ynew=LayoutperDirection.y_new
        layout_ynew.loc[LayoutperDirection.y_new.isnull()]=LayoutperDirection.y0.loc[LayoutperDirection.y_new.isnull()]
        fi.reinitialize(layout_x=np.array(LayoutperDirection.x_new).tolist(),
                        layout_y=np.array(LayoutperDirection.y_new).tolist()) #layout fillna
        AEP_new=fi.get_farm_AEP([freq[ww]],freq_warn=0)
        if plot==1:
            hor_plane_new = fi.calculate_horizontal_plane(height=fi.floris.farm.turbine_definitions[0]['hub_height'],x_bounds=[-3000,3000],y_bounds=[-3000,3000])
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
            ax1.set_xlim(-1920,1920)
            ax0.set_xlim(-1920,1920)
            ax1.set_ylim(-1920,1920)
            ax0.set_ylim(-1920,1920)
            ax1.set_xticks(np.arange(-1440,2800,1440))
            ax0.set_xticks(np.arange(-1440,2800,1440))
            ax0.set_yticks(np.arange(-1440,2800,1440))
            ax1.set_yticks(np.arange(-1440,2800,1440))
            ax0.set_xticklabels(np.arange(-1440/240,2800/240,1440/240))
            ax1.set_xticklabels(np.arange(-1440/240,2800/240,1440/240))
            ax0.set_yticklabels(np.arange(-1440/240,2800/240,1440/240))
            visualize_cut_plane(hor_plane0, ax=ax0)
            im=visualize_cut_plane(hor_plane_new, ax=ax1)
            
            ax0.plot(layout[0],layout[1],'ok')
            ax1.plot(layout[0],layout[1],'ok')
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
            ["OWFL"],
            loc="lower center",
            bbox_to_anchor=(0.9, 0.92),
            ncol=1,
            fontsize=11)
            plt.subplots_adjust(wspace=0.1)
            # plt.savefig(path+'/Figures_WES'+ '/Target_Finallayout_'+'Wind_dir_%s.png'%(wd[ww]),dpi=200)
            # plt.savefig(r'/home/youssef/Documents/Conferences/WESC/Figures/Target_vs_Finallayout_'+'Wind_dir_%s.png'%(wd[ww]),dpi=200)
    delta_E=pd.DataFrame(delta_AEPtotal,columns=['delta_AEP'])
    delta_E['AEP0']=total_AEP0
    delta_E['AEP_new']=total_AEP_new
    delta_E['Wind_dir']=wd
    delta_E['Wind_dir']=90-delta_E['Wind_dir']+180 #changing wind from floris to global coordinates
    delta_E['Wind_dir']%= 360
    delta_E['Wind_dir']=myround(delta_E['Wind_dir'], base=5)
    delta_E.sort_values(by=['delta_AEP'],ignore_index=True,inplace=True)
    return delta_E

def turbines_velocity(wd,ws,freq,layout0,newlayoutdict,delta=1):
    turbines_V0 = pd.DataFrame()
    turbines_V = pd.DataFrame()
    layout=layout0   
    NewLayoutX_delta= newlayoutdict["X"].copy()
    NewLayoutY_delta= newlayoutdict["Y"].copy()
    NewLayoutX_delta=NewLayoutX_delta.reindex(sorted(NewLayoutX_delta.columns), axis=1)
    NewLayoutY_delta=NewLayoutY_delta.reindex(sorted(NewLayoutY_delta.columns), axis=1)
    
    col=pd.DataFrame(NewLayoutX_delta.columns.astype(str).astype(float),columns=['MoorPy_ang'])
    col['Floris_ang'] = 90-col['MoorPy_ang']+180
    col['Floris_ang'] %= 360
    
    NewLayoutX_delta=NewLayoutX_delta.set_axis(col['Floris_ang'],axis=1)
    NewLayoutY_delta=NewLayoutY_delta.set_axis(col['Floris_ang'],axis=1)
    
    NewLayoutX_delta=NewLayoutX_delta.reindex(sorted(NewLayoutX_delta.columns), axis=1)
    NewLayoutY_delta=NewLayoutY_delta.reindex(sorted(NewLayoutY_delta.columns), axis=1)
    NewLayoutX_delta=NewLayoutX_delta.set_axis(wd,axis=1)
    NewLayoutY_delta=NewLayoutY_delta.set_axis(wd,axis=1)
    # pdb.set_trace() 
    
    for ww in range(len(wd)):
        LayoutperDirection=pd.DataFrame(layout).transpose()
        LayoutperDirection.columns=['x0','y0']
        LayoutperDirection['x_new']= LayoutperDirection.x0
        LayoutperDirection['y_new']= LayoutperDirection.y0
        if delta==1:
            LayoutperDirection['x_new']= NewLayoutX_delta.loc[:,wd[ww]]+LayoutperDirection.x0
            LayoutperDirection['y_new']= NewLayoutY_delta.loc[:,wd[ww]]+LayoutperDirection.y0 
        if delta==0:
            # pdb.set_trace() 
            LayoutperDirection['x_new']= NewLayoutX_delta.loc[:,wd[ww]]
            LayoutperDirection['y_new']= NewLayoutY_delta.loc[:,wd[ww]]

        fi.reinitialize(layout_x=np.array(LayoutperDirection.x0).tolist(),
                        layout_y=np.array(LayoutperDirection.y0).tolist(),
                        wind_directions=[wd[ww]],
                        wind_speeds=ws) 

        fi.calculate_wake()
        V=fi.get_turbine_average_velocities()
        V=V.reshape(1,no_turbines)
        V=V.transpose()
        # pdb.set_trace() 
        turbines_V0[wd[ww]] = pd.DataFrame(V)
    
        fi.reinitialize(layout_x=np.array(LayoutperDirection.x_new).tolist(),
                        layout_y=np.array(LayoutperDirection.y_new).tolist(),
                        wind_directions=[wd[ww]],
                        wind_speeds=ws) #layout fillna
        
        fi.calculate_wake()
        V=fi.get_turbine_average_velocities()
        V=V.reshape(1,no_turbines)
        V=V.transpose()
        turbines_V[wd[ww]] = pd.DataFrame(V)
    
    turbines_V0.columns = 270 - turbines_V0.columns
    turbines_V0.columns %= 360
    turbines_V0 = turbines_V0.reindex(sorted(turbines_V0.columns), axis=1)
    turbines_V.columns = 270 - turbines_V.columns
    turbines_V.columns %= 360
    turbines_V = turbines_V.reindex(sorted(turbines_V.columns), axis=1)
        
    return turbines_V0, turbines_V



########### Inputs start  
mindist = 4
spc = 6 # maximum spacing at the beginning of the simulation for boundary creation
shape ='circle'
no_circles_squares =1 
no_turbines = 7 
file_number = f'{mindist}{mindist}{no_turbines}{no_turbines}'
windrose = 'iea'
yawang = 5
Ti = 0.06
angle=1

############################ 
path0=fr"../Results_{windrose}"
path =fr"../Results_{windrose}/Results_{mindist}_{spc}"

opt_baseline_file ='/SNOPTlayoutnew_'+shape+'_'+str(no_turbines)+'T_iea.csv'
opt_dyn_layout_file='/OptimumDynamicfarmlayout_'+shape+'_'+str(no_turbines)+'T_iea_TI_'+str(Ti)+'.csv'
opt_dyn_layout_file_xaxis='/OptimumDynamicfarmlayout_xaxis_'+shape+'_'+str(no_turbines)+'T_iea_TI_'+str(Ti)+'.csv'
opt_dyn_layout_file_yaxis='/OptimumDynamicfarmlayout_yaxis_'+shape+'_'+str(no_turbines)+'T_iea_TI_'+str(Ti)+'.csv'
opt_aep_sum_file ='/AEP_summary_per_winddir_'+shape+'_'+str(no_turbines)+'T_iea_TI_'+str(Ti)+'.csv'
outputfile='/TBMDfinal_Yaw'+ str(yawang) + 'deg'+f'_{shape}_{no_turbines}T_iea.csv'

# Mooring Database inputs
MoorPy_disp= pd.read_csv(path0+'/MooringDatabase_Yaw5deg/WindSpeed_10'+ '/perp_diff_Sway.csv').round(decimals=2) # The perpendicular displacement of each mooring design
MooringPermMatrix= pd.read_csv(path0+'/MooringDatabase_Yaw5deg/WindSpeed_10'+ '/Final_PermMatrix.csv') # the design parameters vof each mooring system
Xdiff_file=path0+'/MooringDatabase_Yaw5deg/WindSpeed_10'+ '/Xdiff_Surge.csv' # The positions achieved by each mooring design in each wind dir (x-axis)
Ydiff_file=path0+'/MooringDatabase_Yaw5deg/WindSpeed_10'+ '/Xdiff_Sway.csv' # The positions achieved by each mooring design in each wind dir (y-axis)

Floris_perp= pd.read_csv(path+ '/SingleWIndDir_Opt'+opt_dyn_layout_file)
Floris_perp.columns=Floris_perp.columns.astype(str).astype(float)
Floris_perp_xaxis= pd.read_csv(path+ '/SingleWIndDir_Opt'+opt_dyn_layout_file_xaxis)
Floris_perp_xaxis.columns=Floris_perp_xaxis.columns.astype(str).astype(float)
Floris_perp_yaxis= pd.read_csv(path+ '/SingleWIndDir_Opt'+opt_dyn_layout_file_yaxis)
Floris_perp_yaxis.columns=Floris_perp_yaxis.columns.astype(str).astype(float)
AEP_angle= pd.read_csv(path+ '/SingleWIndDir_Opt'+opt_aep_sum_file) # I think maybe we dont need to read this file at all and get the FLORIS wind dir from the windrose data only
# pdb.set_trace() 

Xdiff=pd.read_csv(Xdiff_file).round(decimals=2)
Ydiff=pd.read_csv(Ydiff_file).round(decimals=2)

layout_final=pd.read_csv(path+ '/FinalMooringWFDesign' + outputfile,index_col=False)

# Windrose parameters
# Windrose parameters
if windrose =='iea':
    print('############################ iea windrose ############################')
    wd = np.array([0., 22.5, 45., 67.5, 90., 112.5, 135., 157.5, 180., 202.5, 225., 247.5, 270., 292.5, 315., 337.5], dtype=float)
    ws = np.array([10],dtype=float)
    ws_d = np.array([10, 10,  10,  10,    10,   10,    10,   10,    10,    10,    10,   10,   10,    10,   10,    10], dtype=float)
    freq_d = [.025,  .024,  .029,  .036,.063,  .065,  .100,  .122,.063,  .038,  .039,  .083, .213,  .046,  .032,  .022]
    ########################################################
if windrose =='alpha':
    print(' ########################### alpha ventus windrose ############################')
    wd = np.array([0., 22.5, 45., 67.5, 90., 112.5, 135., 157.5, 180., 202.5, 225., 247.5, 270., 292.5, 315., 337.5], dtype=float)
    ws = np.array([10],dtype=float)
    ws_d = np.array([10, 10,  10,  10,    10,   10,    10,   10,    10,    10,    10,   10,   10,    10,   10,    10], dtype=float)
    freq_d = [0.0313,0.0402,0.0375, 0.0568,0.0558,0.0608, 0.0424,0.0564,0.0555, 0.1114,0.0932,0.1114, 0.0722,
          0.0743,0.0500, 0.0508]
    # #######################################################

freq=np.transpose(np.array([freq_d]))

# N_circles=1
# D=240
# spc=5
# layout_base=ConcentricCirclesLayout(N_circles,spc,D)
D=240
############## Circle  ############## 
if shape=='circle':
    N_circles=no_circles_squares
    spacing=spc*D
    angles=np.arange(0,360+angle,angle)
    boundaries_x=np.round(((N_circles*spacing)+0.01*D)*np.cos(np.radians(angles)))
    boundaries_y=np.round(((N_circles*spacing)+0.01*D)*np.sin(np.radians(angles)))
    boundaries = [[x,y] for x, y in zip(boundaries_x, boundaries_y)]
    # boundaries=sort_boundaries(boundaries)
    layout_base=ConcentricCirclesLayout(N_circles,spc,D)    


############## Square  ############## 
## D=240
if shape=='square':
    N_squares=no_circles_squares
    spc=spc*D
    y_coordinates=np.tile(np.arange(-spc*N_squares,spc*N_squares+spc,spc),1+N_squares*2)
    x_coordinates=np.repeat(np.arange(-spc*N_squares,spc*N_squares+spc,spc),1+N_squares*2)
    x_coordinates=np.round(x_coordinates)
    y_coordinates=np.round(y_coordinates)
    layout_base= (x_coordinates.tolist(), y_coordinates.tolist())
    





# wind farm layout same as the one used in Pyoptsparsetry1
layoutdf=pd.read_csv(path+ '/BaselineOptimization' +opt_baseline_file)

layout0=(np.array(layoutdf.x).tolist(),
        np.array(layoutdf.y).tolist())
# wind tuebine parameters and wake model
fi = FlorisInterface('../Inputs/gch15MW.yaml')
fi.floris.flow_field.turbulence_intensity=Ti

pdb.set_trace() 

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
# plt.savefig(path+  '/Figures'+ '/Baseline_25T.png',dpi=200)
plt.show()





################ Optimized layout ###############
mydpi=96
fig= plt.figure(figsize=(5,5),dpi=mydpi)
fontsize = 13
ax = fig.add_subplot()
ax.plot(np.array(layout_base[0])/240, np.array(layout_base[1])/240, "ob")
ax.plot(np.array(layout0[0])/240, np.array(layout0[1])/240, "xr")
ax.set_xlabel("x/D", fontsize=fontsize)
ax.set_ylabel("y/D", fontsize=fontsize,loc='center')
# plt.axis("equal")
# ax.set_yticks([-10,-5,0,5,10])
# ax.set_xticks([-7,-5,0,5,7])
# boundaries_x=np.round((min(layout_base[0])-0.01*D,min(layout_base[0])-0.01*D,
#                         max(layout_base[0])+0.01*D,max(layout_base[0])+0.01*D,
#                         min(layout_base[0])-0.01*D),decimals=2)
# boundaries_y=np.round((min(layout_base[0])-0.01*D,max(layout_base[0])+0.01*D,
#                         max(layout_base[0])+0.01*D,min(layout_base[0])-0.01*D,
#                         min(layout_base[0])-0.01*D),decimals=2)
ax.xaxis.set_tick_params( labelsize=fontsize)
ax.yaxis.set_tick_params( labelsize=fontsize)
# ax.plot(boundaries_x/240, boundaries_y/240, "-k")
ax.legend(
    ["Baseline", "Optimized"],
    loc="lower center",
    bbox_to_anchor=(0.5, 1.01),
    ncol=2,
    fontsize=fontsize,
        )
ax.yaxis.set_label_coords(-0.1, 0.5)
plt.grid()
# plt.savefig(path+'/Figures_WES/Optimized_9T_iea_wobound.png',dpi=200)
plt.show()


################ Windrose constant wind speed ###############
mydpi=96
fontsize = 13
legend_kwargs=dict(ncol=1, loc='best', fancybox=True,bbox_to_anchor=(0.4, 0.55, 0.5, 0.5))
fig = plt.figure(figsize=(650/mydpi,450/mydpi),dpi=mydpi)
ax = fig.add_subplot(projection='polar')
# ax.set_yticklabels([])
ax.xaxis.set_tick_params( labelsize=fontsize)
ax.yaxis.set_tick_params( labelsize=fontsize)
ax.set_rticks([0.05,0.1,0.15,0.2,0.25])
ax.set_rlabel_position(65.5)
Rose_onews = WindRose()
Data_Rose_onews = Rose_onews.make_wind_rose_from_user_dist(np.array(wd), np.array(ws_d), np.array(freq_d), wd=np.array(wd),ws=np.arange(3,26,1))
Rose_onews.plot_wind_rose(ws_right_edges=np.array([5, 10, 15, 25]),ax=ax,wd_bins=np.array(wd),legend_kwargs=legend_kwargs,color_map='hot')
ax.legend_.remove()
# plt.legend(['10 m/s'],bbox_to_anchor=(0.4, 0.55, 0.5, 0.5))
# plt.savefig(path+'/Figures_WES/windroseconstwind.png',dpi=200)
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
Windrose_real=make_wind_rose_from_weibull(wd, freq_d, wd=wd, ws=np.arange(3,26,1))
Rose = WindRose()
Data_Rose= Rose.make_wind_rose_from_user_dist(np.array(Windrose_real.wd), np.array(Windrose_real.ws),
                                                    np.array( Windrose_real.freq_val), wd=np.array(wd),ws=np.arange(3,26,1))
Rose.plot_wind_rose(ws_right_edges=np.array([10,20]),ax=ax,wd_bins=np.array(wd),legend_kwargs=legend_kwargs)
plt.legend(['> 10 m/s','< 10 m/s'],bbox_to_anchor=(0.4, 0.55, 0.5, 0.5))
# plt.savefig(r'D:\Moordyn\moorpy\MoorPy\Windenergy_journal_fig\windrose.png',dpi=200)
plt.show()


############### 

fi.reinitialize(layout=layout0,
                wind_directions=wd,
                wind_speeds=ws) 
# fi.reinitialize(layout=layout0)
AEP_initial=fi.get_farm_AEP(freq)
AEP_initial_nowake=fi.get_farm_AEP(freq,no_wake=True)
wake_losses=(AEP_initial-AEP_initial_nowake)/AEP_initial_nowake*100
print(wake_losses)
fi.calculate_wake()
# V_baseline=fi.get_turbine_average_velocities()
# V_baseline=V_baseline.reshape(16,9)
# V_baseline=V_baseline.transpose()
# wd_FAST=270-wd
# wd_FAST%= 360
# V_baseline=pd.DataFrame(V_baseline, columns=wd_FAST)
# V_baseline.reindex(sorted(V_baseline.columns), axis=1)
# fi.reinitialize(layout=layout_base)
# AEP_base=fi.get_farm_AEP(freq)

############### Targeted layout ###############
TargetLayoutdict={"X":Floris_perp_xaxis,"Y":Floris_perp_yaxis}
delta_E_target = delta_energy(wd,ws,freq,layout0,TargetLayoutdict,0,AEP_initial,delta=0)    
delta_E_target.sort_values(by=['Wind_dir'],inplace=True, ignore_index=True)
V0,Vt=turbines_velocity(wd,ws,freq,layout0,TargetLayoutdict,delta=0)
gain0=sum(delta_E_target['delta_AEP'])
print('gain0=')
print(gain0) 

# V0.to_csv(path+ '/FinalMooringWFDesign' +'/turbinevelocities_baseline.csv', index=False)
# Vt.to_csv(path+ '/FinalMooringWFDesign' + '/turbinevelocities_target.csv', index=False)

# pdb.set_trace() 
############### Final layout ###############
Layoutdict=layoutfromMoorSystem(layout_final,Floris_perp,Xdiff,Ydiff)
delta_E_final = delta_energy(wd,ws,freq,layout0,Layoutdict,0,AEP_initial,delta=1) 
delta_E_final.sort_values(by=['Wind_dir'],inplace=True, ignore_index=True)   
V0,Vf=turbines_velocity(wd,ws,freq,layout0,Layoutdict,delta=1) 
gain0=sum(delta_E_final['delta_AEP'])
print('gain0=')
print(gain0)       

# Vf.to_csv(path+ '/FinalMooringWFDesign/Yaw' + str(yawang) + 'deg' +'/turbinevelocities_final.csv', index=False)

delta_E_finalplot=delta_E_final.copy()
delta_E_finalplot['Wind_dir']=90-delta_E_finalplot['Wind_dir']+180+360
delta_E_finalplot['Wind_dir']%=360
delta_E_finalplot.sort_values(by=['Wind_dir'],inplace=True, ignore_index=True)  
delta_E_finalplot['Wind_dir']=wd
delta_E_finalplot.to_csv(path+ '/FinalMooringWFDesign' + '/delta_E_finalplot.csv', index=False)

delta_E_targetplot=delta_E_target.copy()
delta_E_targetplot['Wind_dir']=90-delta_E_targetplot['Wind_dir']+180+360
delta_E_targetplot['Wind_dir']%=360
delta_E_targetplot.sort_values(by=['Wind_dir'],inplace=True, ignore_index=True)  
delta_E_targetplot['Wind_dir']=wd
delta_E_targetplot.to_csv(path+ '/FinalMooringWFDesign' + '/delta_E_targetplot.csv', index=False)


############### Final gain ###############
mydpi=96
fig = plt.figure(figsize=(650/mydpi,450/mydpi),dpi=mydpi)
ax = fig.add_subplot(projection='polar')
ax.bar(np.deg2rad(delta_E_finalplot['Wind_dir']),delta_E_finalplot['delta_AEP'],width= np.radians(22.5),
       edgecolor="k",color='g')
ax.set_theta_direction(-1)
ax.set_theta_offset(np.pi / 2.0)
# ax.set_rticks([-0.2,0.2,0.4,0.6])
fontsize=12
ax.xaxis.set_tick_params( labelsize=fontsize)
ax.yaxis.set_tick_params( labelsize=fontsize)
ax.set_theta_zero_location("N")
ax.set_xticklabels(["N", "NE", "E", "SE", "S", "SW", "W", "NW"])
# ax.set_yticklabels(["-0.2","0.2", "0.4", "0.6%"])
ax.set_title('Final energy gain',fontsize=14)
ax.set_rlabel_position(65.5)
# plt.savefig(r'/home/youssef/Documents/Conferences/WESC/Figures/Final_E_gain.png',dpi=200)
plt.show()


############### Target gain ###############
fig = plt.figure(figsize=(650/mydpi,450/mydpi),dpi=mydpi)
ax = fig.add_subplot(projection='polar')
ax.bar(np.deg2rad(delta_E_targetplot['Wind_dir']),delta_E_targetplot['delta_AEP'],width= np.radians(22.5),
       edgecolor="k",color='g')
ax.set_theta_direction(-1)
ax.set_theta_offset(np.pi / 2.0)
ax.set_rticks([0.4,0.8,1.2,1.6])
ax.set_theta_zero_location("N")
ax.xaxis.set_tick_params( labelsize=fontsize)
ax.yaxis.set_tick_params( labelsize=fontsize)
ax.set_xticklabels(["N", "NE", "E", "SE", "S", "SW", "W", "NW"])
ax.set_yticklabels([ "0.4","0.8","1.2","1.6%"])
ax.set_title('Targeted energy gain',fontsize=14)
ax.set_rlabel_position(65.5)
# plt.savefig(r'/home/youssef/Documents/Conferences/WESC/Figures/Targeted_E_gain.png',dpi=200)
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
ax0.set_rticks([0.25,0.5,1])
fontsize=12
ax0.xaxis.set_tick_params( labelsize=fontsize)
ax0.yaxis.set_tick_params( labelsize=fontsize)
ax0.set_theta_zero_location("N")
ax0.set_xticklabels(["N", "NE", "E", "SE", "S", "SW", "W", "NW"])
ax0.set_yticklabels(["0.25","0.5", "1%"])
ax0.set_title('Final energy gain',fontsize=14)
ax0.set_rlabel_position(65.5)
ax1.bar(np.deg2rad(delta_E_targetplot['Wind_dir']),delta_E_targetplot['delta_AEP'],width= np.radians(22.5),
       edgecolor="k",color='g')
ax1.set_theta_direction(-1)
ax1.set_theta_offset(np.pi / 2.0)
ax1.set_rticks([0.25,0.5,1])
ax1.set_theta_zero_location("N")
ax1.xaxis.set_tick_params( labelsize=fontsize)
ax1.yaxis.set_tick_params( labelsize=fontsize)
ax1.set_xticklabels(["N", "NE", "E", "SE", "S", "SW", "W", "NW"])
ax1.set_yticklabels(["0.25","0.5", "1%"])
ax1.set_title('Targeted energy gain',fontsize=14)
ax1.set_rlabel_position(65.5)
# plt.savefig(r'D:\Moordyn\moorpy\MoorPy\Windenergy_journal_fig\Targeted_Final_E_gain.png',dpi=200)
plt.show()


################ All Baseline layouts ###############
# mydpi=96
# fig = plt.figure(figsize=(9,9),dpi=mydpi)
# ax0 = fig.add_subplot(221)
# ax1 = fig.add_subplot(222)
# ax2 = fig.add_subplot(223)
# ax3 = fig.add_subplot(224)

# N_circles=1
# D=240
# spc=5
# spacing=5*D
# spc=5
# angles=np.arange(0,360,1)
# layout_base=ConcentricCirclesLayout(N_circles,spc,D)
# boundaries_x=np.round(((N_circles*spacing)+0.01*D)*np.cos(np.radians(angles)))
# boundaries_y=np.round(((N_circles*spacing)+0.01*D)*np.sin(np.radians(angles)))
# fontsize=12
# ax0.plot(np.array(layout_base[0])/240, np.array(layout_base[1])/240, "ob")
# ax0.plot(boundaries_x/240, boundaries_y/240, "-k")
# ax0.set_yticks([-15,-10,-5,0,5,10,15])
# ax0.set_xticks([-15,-10,-5,0,5,10,15])
# ax0.xaxis.set_tick_params( labelsize=fontsize)
# ax0.yaxis.set_tick_params( labelsize=fontsize)
# ax0.set_ylabel("y/D", fontsize=fontsize,loc='center')
# ax0.grid()
# ax0.legend(
#     ["Turbines", "Boundaries"],
#     ncol=1,
#     fontsize=fontsize,
#         )
# ax0.set_title('7 turbines',fontsize=14)


# N_circles=3
# D=240
# spc=5
# spacing=5*D
# spc=5
# angles=np.arange(0,360,1)
# layout_base=ConcentricCirclesLayout(N_circles,spc,D)
# boundaries_x=np.round(((N_circles*spacing)+0.01*D)*np.cos(np.radians(angles)))
# boundaries_y=np.round(((N_circles*spacing)+0.01*D)*np.sin(np.radians(angles)))
# fontsize=12
# ax1.plot(np.array(layout_base[0])/240, np.array(layout_base[1])/240, "ob")
# ax1.plot(boundaries_x/240, boundaries_y/240, "-k")
# ax1.set_yticks([-15,-10,-5,0,5,10,15])
# ax1.set_xticks([-15,-10,-5,0,5,10,15])
# ax1.xaxis.set_tick_params( labelsize=fontsize)
# ax1.yaxis.set_tick_params( labelsize=fontsize)
# ax1.grid()
# ax1.set_title('37 turbines',fontsize=14)

# N_squares=1
# D=240
# # spacing=5*D
# spc=5
# spc=spc*D
# y_coordinates=np.tile(np.arange(-spc*N_squares,spc*N_squares+spc,spc),1+N_squares*2)
# x_coordinates=np.repeat(np.arange(-spc*N_squares,spc*N_squares+spc,spc),1+N_squares*2)
# x_coordinates=np.round(x_coordinates)
# y_coordinates=np.round(y_coordinates)
# layout_base= (x_coordinates.tolist(), y_coordinates.tolist())
# boundaries_x=np.round((min(layout_base[0])-0.01*D,min(layout_base[0])-0.01*D,
#                         max(layout_base[0])+0.01*D,max(layout_base[0])+0.01*D,
#                         min(layout_base[0])-0.01*D),decimals=2)
# boundaries_y=np.round((min(layout_base[0])-0.01*D,max(layout_base[0])+0.01*D,
#                         max(layout_base[0])+0.01*D,min(layout_base[0])-0.01*D,
#                         min(layout_base[0])-0.01*D),decimals=2)
# ax2.plot(np.array(layout_base[0])/240, np.array(layout_base[1])/240, "ob")
# ax2.plot(boundaries_x/240, boundaries_y/240, "-k")
# ax2.set_yticks([-15,-10,-5,0,5,10,15])
# ax2.set_xticks([-15,-10,-5,0,5,10,15])
# ax2.xaxis.set_tick_params( labelsize=fontsize)
# ax2.yaxis.set_tick_params( labelsize=fontsize)
# ax2.set_xlabel("x/D", fontsize=fontsize)
# ax2.set_ylabel("y/D", fontsize=fontsize,loc='center')
# ax2.grid()
# ax2.set_title('9 turbines',fontsize=14)

# N_squares=2
# D=240
# # spacing=5*D
# spc=5
# spc=spc*D
# y_coordinates=np.tile(np.arange(-spc*N_squares,spc*N_squares+spc,spc),1+N_squares*2)
# x_coordinates=np.repeat(np.arange(-spc*N_squares,spc*N_squares+spc,spc),1+N_squares*2)
# x_coordinates=np.round(x_coordinates)
# y_coordinates=np.round(y_coordinates)
# layout_base= (x_coordinates.tolist(), y_coordinates.tolist())
# boundaries_x=np.round((min(layout_base[0])-0.01*D,min(layout_base[0])-0.01*D,
#                         max(layout_base[0])+0.01*D,max(layout_base[0])+0.01*D,
#                         min(layout_base[0])-0.01*D),decimals=2)
# boundaries_y=np.round((min(layout_base[0])-0.01*D,max(layout_base[0])+0.01*D,
#                         max(layout_base[0])+0.01*D,min(layout_base[0])-0.01*D,
#                         min(layout_base[0])-0.01*D),decimals=2)
# ax3.plot(np.array(layout_base[0])/240, np.array(layout_base[1])/240, "ob")
# ax3.plot(boundaries_x/240, boundaries_y/240, "-k")
# ax3.set_yticks([-15,-10,-5,0,5,10,15])
# ax3.set_xticks([-15,-10,-5,0,5,10,15])
# ax3.xaxis.set_tick_params( labelsize=fontsize)
# ax3.yaxis.set_tick_params( labelsize=fontsize)
# ax3.set_xlabel("x/D", fontsize=fontsize)
# ax3.grid()
# ax3.set_title('25 turbines',fontsize=14)

# plt.savefig(path+  '/Figures'+ '/Baseline_alllayouts_iea.png',dpi=200)
# plt.show()



################ All Optimized layouts ###############
# mydpi=96
# fig = plt.figure(figsize=(9,9),dpi=mydpi)
# ax0 = fig.add_subplot(221)
# ax1 = fig.add_subplot(222)
# ax2 = fig.add_subplot(223)
# ax3 = fig.add_subplot(224)

# N_circles=1
# D=240
# spc=5
# spacing=5*D
# spc=5
# angles=np.arange(0,360,1)
# layout_base=ConcentricCirclesLayout(N_circles,spc,D)
# boundaries_x=np.round(((N_circles*spacing)+0.01*D)*np.cos(np.radians(angles)))
# boundaries_y=np.round(((N_circles*spacing)+0.01*D)*np.sin(np.radians(angles)))
# layoutdf=pd.read_csv(path+ '/BaselineOptimization' +'/SNOPTlayoutnew_circle_7T_iea_10degrees.csv')
# layout0=(np.array(layoutdf.x).tolist(),
#         np.array(layoutdf.y).tolist())
# fontsize=12
# ax0.plot(np.array(layout_base[0])/240, np.array(layout_base[1])/240, "ob")
# ax0.plot(np.array(layout0[0])/240, np.array(layout0[1])/240, "xr")
# ax0.plot(boundaries_x/240, boundaries_y/240, "-k")
# ax0.set_yticks([-15,-10,-5,0,5,10,15])
# ax0.set_xticks([-15,-10,-5,0,5,10,15])
# ax0.xaxis.set_tick_params( labelsize=fontsize)
# ax0.yaxis.set_tick_params( labelsize=fontsize)
# ax0.set_ylabel("y/D", fontsize=fontsize,loc='center')
# ax0.grid()
# ax0.legend(
#     ["Baseline", "Optimized","Boundaries"],
#     ncol=1,
#     fontsize=fontsize,
#         )
# ax0.set_title('7 turbines',fontsize=14)


# N_circles=3
# D=240
# spc=5
# spacing=5*D
# spc=5
# angles=np.arange(0,360,1)
# layout_base=ConcentricCirclesLayout(N_circles,spc,D)
# boundaries_x=np.round(((N_circles*spacing)+0.01*D)*np.cos(np.radians(angles)))
# boundaries_y=np.round(((N_circles*spacing)+0.01*D)*np.sin(np.radians(angles)))
# layoutdf=pd.read_csv(path+ '/BaselineOptimization' +'/SNOPTlayoutnew_circle_37T_iea37.csv')
# layout0=(np.array(layoutdf.x).tolist(),
#         np.array(layoutdf.y).tolist())
# ax1.plot(np.array(layout_base[0])/240, np.array(layout_base[1])/240, "ob")
# ax1.plot(np.array(layout0[0])/240, np.array(layout0[1])/240, "xr")
# ax1.plot(boundaries_x/240, boundaries_y/240, "-k")
# ax1.set_yticks([-15,-10,-5,0,5,10,15])
# ax1.set_xticks([-15,-10,-5,0,5,10,15])
# ax1.xaxis.set_tick_params( labelsize=fontsize)
# ax1.yaxis.set_tick_params( labelsize=fontsize)
# ax1.grid()
# ax1.set_title('37 turbines',fontsize=14)

# N_squares=1
# D=240
# # spacing=5*D
# spc=5
# spc=spc*D
# y_coordinates=np.tile(np.arange(-spc*N_squares,spc*N_squares+spc,spc),1+N_squares*2)
# x_coordinates=np.repeat(np.arange(-spc*N_squares,spc*N_squares+spc,spc),1+N_squares*2)
# x_coordinates=np.round(x_coordinates)
# y_coordinates=np.round(y_coordinates)
# layout_base= (x_coordinates.tolist(), y_coordinates.tolist())
# boundaries_x=np.round((min(layout_base[0])-0.01*D,min(layout_base[0])-0.01*D,
#                         max(layout_base[0])+0.01*D,max(layout_base[0])+0.01*D,
#                         min(layout_base[0])-0.01*D),decimals=2)
# boundaries_y=np.round((min(layout_base[0])-0.01*D,max(layout_base[0])+0.01*D,
#                         max(layout_base[0])+0.01*D,min(layout_base[0])-0.01*D,
#                         min(layout_base[0])-0.01*D),decimals=2)
# layoutdf=pd.read_csv(path+ '/BaselineOptimization' +'/SNOPTlayoutnew_square_9T_iea37.csv')
# layout0=(np.array(layoutdf.x).tolist(),
#         np.array(layoutdf.y).tolist())
# ax2.plot(np.array(layout_base[0])/240, np.array(layout_base[1])/240, "ob")
# ax2.plot(np.array(layout0[0])/240, np.array(layout0[1])/240, "rx")
# ax2.plot(boundaries_x/240, boundaries_y/240, "-k")
# ax2.set_yticks([-15,-10,-5,0,5,10,15])
# ax2.set_xticks([-15,-10,-5,0,5,10,15])
# ax2.xaxis.set_tick_params( labelsize=fontsize)
# ax2.yaxis.set_tick_params( labelsize=fontsize)
# ax2.set_xlabel("x/D", fontsize=fontsize)
# ax2.set_ylabel("y/D", fontsize=fontsize,loc='center')
# ax2.grid()
# ax2.set_title('9 turbines',fontsize=14)

# N_squares=2
# D=240
# # spacing=5*D
# spc=5
# spc=spc*D
# y_coordinates=np.tile(np.arange(-spc*N_squares,spc*N_squares+spc,spc),1+N_squares*2)
# x_coordinates=np.repeat(np.arange(-spc*N_squares,spc*N_squares+spc,spc),1+N_squares*2)
# x_coordinates=np.round(x_coordinates)
# y_coordinates=np.round(y_coordinates)
# layout_base= (x_coordinates.tolist(), y_coordinates.tolist())
# boundaries_x=np.round((min(layout_base[0])-0.01*D,min(layout_base[0])-0.01*D,
#                         max(layout_base[0])+0.01*D,max(layout_base[0])+0.01*D,
#                         min(layout_base[0])-0.01*D),decimals=2)
# boundaries_y=np.round((min(layout_base[0])-0.01*D,max(layout_base[0])+0.01*D,
#                         max(layout_base[0])+0.01*D,min(layout_base[0])-0.01*D,
#                         min(layout_base[0])-0.01*D),decimals=2)
# layoutdf=pd.read_csv(path+ '/BaselineOptimization' +'/SNOPTlayoutnew_square_25T_iea37.csv')
# layout0=(np.array(layoutdf.x).tolist(),
#         np.array(layoutdf.y).tolist())
# ax3.plot(np.array(layout_base[0])/240, np.array(layout_base[1])/240, "ob")
# ax3.plot(np.array(layout0[0])/240, np.array(layout0[1])/240, "rx")
# ax3.plot(boundaries_x/240, boundaries_y/240, "-k")
# ax3.set_yticks([-15,-10,-5,0,5,10,15])
# ax3.set_xticks([-15,-10,-5,0,5,10,15])
# ax3.xaxis.set_tick_params( labelsize=fontsize)
# ax3.yaxis.set_tick_params( labelsize=fontsize)
# ax3.set_xlabel("x/D", fontsize=fontsize)
# ax3.grid()
# ax3.set_title('25 turbines',fontsize=14)

# plt.savefig(path+  '/Figures'+ '/Optimized_alllayouts_iea.png',dpi=200)
# plt.show()





delta_energyplot(wd,ws,freq,TargetLayoutdict,Layoutdict,1,AEP_initial,layout0,delta=1)


mydpi=96
fig = plt.figure(figsize=(4.5,5),dpi=mydpi)
ax = fig.add_subplot(projection='polar')
aero_angle1=np.arange(0,360,5)+55
# ax.plot(np.deg2rad(aero_angle1)
#         ,np.abs(MoorPy_disp.loc[3000,:]),'.b')

aero_angle=Floris_perp.columns.values
aero_angle %= 360
columnsname=Floris_perp.columns.values.astype(int)-55
columnsname %= 360
# ax.plot(np.deg2rad(aero_angle)
#         ,np.abs(MoorPy_disp.loc[3000,columnsname.astype(str)]),'.g')


# rnd=np.random.randint(110,size=(18))
# ax.plot(np.append(np.arange(0,2*np.pi,np.pi/9),0)
#         ,np.append(rnd,rnd[0]),'r')
# ax.plot(np.deg2rad(Floris_perp.columns.values),
#         abs(Floris_perp.loc[0,:]),'.r')
# ax.set_theta_direction(-1)
# ax.set_theta_offset(np.pi / 2.0)

diff=abs(-np.abs(MoorPy_disp.loc[3000,columnsname.astype(str)].values)+abs(Floris_perp.loc[0,:].values))
ax.plot(np.deg2rad(aero_angle)
        ,diff,'.b')
fontsize=12
ax.set_rticks([40,80,120,160])
ax.set_yticklabels([40,80,120,160])
ax.set_theta_zero_location("N")
ax.xaxis.set_tick_params( labelsize=fontsize)
ax.yaxis.set_tick_params( labelsize=fontsize)
ax.set_xticklabels(["N", "NE", "E", "SE", "S", "SW", "W", "NW"])
ax.legend(
    ["diff crosswind displacement"],
    loc="lower center",
    bbox_to_anchor=(0.5, 1.08),
    ncol=2,
    fontsize=fontsize,
        )
# ax.set_title('Targeted energy gain',fontsize=14)
# ax.set_rlabel_position(65.5)
plt.savefig(r'Results_v2/MooringDatabase/FullResults/WindSpeed_10/Figures/diff.png',dpi=200)
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
Rose_onews = wind_rose.WindRose()
Data_Rose_onews = Rose_onews.make_wind_rose_from_user_dist(np.array(wd), np.array(ws), np.array(freq), wd=np.array(wd),ws=np.arange(3,26,1))
Rose_onews.plot_wind_rose(ws_right_edges=np.array([5, 10, 15, 25]),ax=ax1,wd_bins=np.array(wd),legend_kwargs=legend_kwargs)
ax1.legend_.remove()
# plt.savefig(r'D:\Moordyn\moorpy\MoorPy\Windenergy_journal_fig\Baseline_windrose_together.png',dpi=200)
plt.show()


import pandas as pd
import numpy as np
import pdb
from floris.tools import FlorisInterface
from floris.tools import WindRose
from floris.tools.visualization import visualize_cut_plane
import matplotlib.pyplot as plt
import math
import time
fi = FlorisInterface('/home/youssef/Documents/FLORIS/gch15MW.yaml')
# wd=[0., 22.5, 45., 67.5, 90., 112.5, 135., 157.5, 180., 202.5, 225., 247.5, 270., 292.5, 315., 337.5]
wd=[0.]
ws=[10]
# ws_d=[16, 16,  16,  16,    16,   16,    16,   16,    16,    16,    16,   16,   16,    16,   16,    16]
# freq_d= [.025,  .024,  .029,  .036,.063,  .065,  .100,  .122,.063,  .038,  .039,  .083, .213,  .046,  .032,  .022]
freq_d= [.025]
########################################################
freq=np.transpose(np.array([freq_d]))
layout_try=(np.array([0]).tolist(),
                np.array([0]).tolist())
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
t0=time.time()
fi.reinitialize(layout=layout0,
                wind_directions=wd,
                wind_speeds=ws )
AEP_try=fi.get_farm_AEP(freq)*10**-9
t1=time.time()
