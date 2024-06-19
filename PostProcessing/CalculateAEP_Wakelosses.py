#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 02:28:59 2024

@author: youssef
"""



import pandas as pd
import numpy as np
import pdb
from floris.tools import FlorisInterface
from floris.tools import WindRose
from floris.tools.visualization import visualize_cut_plane
import matplotlib.pyplot as plt
from scipy.interpolate import NearestNDInterpolator
import os
import seaborn as sns

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size" :10,
    "axes.titlesize":10,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8
})

def myround(x, base=5):
    return base * round(x/base)

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

def weibull(x, k = [2.392578, 2.447266, 2.412109, 2.591797, 2.755859, 2.595703,
             2.583984, 2.548828, 2.470703, 2.607422, 2.626953, 2.326172]
, lam = [9.176929, 9.782334, 9.531809, 9.909545, 10.04269, 9.593921,
             9.584007, 10.51499, 11.39895, 11.68746, 11.63732, 10.08803]
):
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
        Weibullall=np.array([])
        for x in range (3,26):
            Weibull = (k / lam) * (x / lam) ** (k - 1) * np.exp(-((x / lam) ** k))
            Weibullall = np.concatenate((Weibullall,Weibull ))
        return Weibullall


def make_wind_rose_from_weibull(winddir,freqperdir, wd=np.arange(0, 360, 5.0), ws=np.arange(0, 26, 1.0),weibull_k= 0.5 , weibull_lam = 0.5):  

    wind_dir =winddir
    freq_dir = freqperdir
# 
    # pdb.set_trace()
    
    freq_wd = np.interp(wd, wind_dir, freq_dir)
    # freq_ws = weibull(ws, k=weibull_k, lam=weibull_lam)

    freq_tot = np.zeros(len(wd) * len(ws))
    wd_tot = np.zeros(len(wd) * len(ws))
    ws_tot = np.zeros(len(wd) * len(ws))

    count = 0
    for i in range(len(wd)):
        for j in range(len(ws)):
            wd_tot[count] = wd[i]
            ws_tot[count] = ws[j]

            freq_tot[count] = freq_wd[i]*((weibull_k[i] / weibull_lam[i]) * (ws[j] / weibull_lam[i]) ** (weibull_k[i] - 1) * np.exp(-((ws[j] / weibull_lam[i]) ** weibull_k[i])))
            count = count + 1
        # pdb.set_trace()

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
    total_AEP_new=[]
    for ww in range(len(wd)):

        LayoutperDirection=pd.DataFrame()
        LayoutperDirection['x_new']= NewLayoutX_delta.loc[:,wd[ww]]+layout[0]
        LayoutperDirection['y_new']= NewLayoutY_delta.loc[:,wd[ww]]+layout[1]  
  
        AEP0 = AEP_dir[str(wd[ww])].copy()

        fi.reinitialize(layout_x=np.array(LayoutperDirection.x_new).tolist(),
                        layout_y=np.array(LayoutperDirection.y_new).tolist(),
                        wind_directions=[wd[ww]],
                        wind_speeds=ws) #layout fillna

        AEP_new=fi.get_farm_AEP([freq[ww]],freq_warn=0)
        delta_AEP=(AEP_new-AEP0)/(AEP_initial)*100
        delta_AEPtotal.append(delta_AEP)
        total_AEP_new.append(AEP_new)

    delta_E=pd.DataFrame(delta_AEPtotal,columns=['delta_AEP'])
    delta_E['Wind_dir']=wd
    delta_E['AEP_new']=total_AEP_new
    delta_E.sort_values(by=['Wind_dir'],ignore_index=True,inplace=True)

    return delta_E

def layoutfromMoorSystem(TurbineMoorDesign,wd,Xdiff,Ydiff): 
    
    wind_dir= np.array(wd)-270
    wind_dir%= 360
    wind_dir= np.round(wind_dir/5)*5
    wind_dir.sort()
    TurbineLayoutX=pd.DataFrame()
    TurbineLayoutY=pd.DataFrame()
    NewLayoutX=pd.DataFrame([],columns=wind_dir)
    NewLayoutY=pd.DataFrame([],columns=wind_dir)
    for i in range(len(Xdiff)):
        # pdb.set_trace()
        TurbineLayoutX.loc[0,Xdiff.columns.values]=Xdiff.iloc[i,:].transpose().values
        TurbineLayoutY.loc[0,Ydiff.columns.values]=Ydiff.iloc[i,:].transpose().values
        col=pd.DataFrame(TurbineLayoutX.columns.astype(str).astype(float),columns=['Floris_ang'])
        col['Floris_ang'] %= 360
        col['MoorPy_ang'] = TurbineLayoutX.columns
        TurbineLayoutX.rename(columns=col.set_index('MoorPy_ang')['Floris_ang'],inplace=True)
        TurbineLayoutY.rename(columns=col.set_index('MoorPy_ang')['Floris_ang'],inplace=True)
        # pdb.set_trace()
        TurbineLayoutX=TurbineLayoutX.reindex(columns=wind_dir)
        TurbineLayoutY=TurbineLayoutY.reindex(columns=wind_dir)
        # pdb.set_trace()
        NewLayoutX.loc[i,:]=TurbineLayoutX.values
        NewLayoutY.loc[i,:]=TurbineLayoutY.values
        TurbineLayoutX=pd.DataFrame()
        TurbineLayoutY=pd.DataFrame()    
    NewLayoutdict={"X":NewLayoutX,"Y":NewLayoutY}
    return NewLayoutdict

def turbine_velocity(wd,ws,freq,layout0,newlayoutdict):
    # t0=time.time()
    turbines_V = pd.DataFrame()
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
    
    for ww in range(len(wd)):
        LayoutperDirection=pd.DataFrame()
        LayoutperDirection['x_new']= NewLayoutX_delta.loc[:,wd[ww]]+layout[0]
        LayoutperDirection['y_new']= NewLayoutY_delta.loc[:,wd[ww]]+layout[1]  


        fi.reinitialize(layout_x=np.array(LayoutperDirection.x_new).tolist(),
                        layout_y=np.array(LayoutperDirection.y_new).tolist(),
                        wind_directions=[wd[ww]],
                        wind_speeds=ws) #layout fillna

        # pdb.set_trace()
        fi.calculate_wake()
        V=fi.get_turbine_average_velocities()
        # pdb.set_trace()
        V=V.reshape(1,len(LayoutperDirection))
        V=V.transpose()
        turbines_V[wd[ww]] = pd.DataFrame(V)

    turbines_V.columns = 270 - turbines_V.columns
    turbines_V.columns %= 360
    turbines_V = turbines_V.reindex(sorted(turbines_V.columns), axis=1)

    return turbines_V

def Floating_layout(wd,layout0,newlayoutdict):
    # t0=time.time()

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
    LayoutperDirection = pd.DataFrame()
    for ww in range(len(wd)):
        
        LayoutperDirection[f'x_new_{wd[ww]}'] = NewLayoutX_delta.loc[:,wd[ww]]+layout[0]
        LayoutperDirection[f'y_new_{wd[ww]}'] = NewLayoutY_delta.loc[:,wd[ww]]+layout[1] 
    return LayoutperDirection
############ PhD input ################ 
Ti = 0.06
tolerance=50
MooringDatabase= 'L26' 
path0=r"../Results_iea"
path =r"../Results_iea/Results_6_6"
opt_dyn_layout_file='/OptimumDynamicfarmlayout_square_9T_iea_TI_0.06.csv'
opt_aep_sum_file ='/AEP_summary_per_winddir_square_9T_iea_TI_0.06.csv'

outputfile = f'/TBMDfinal_tol_{tolerance}_filter_designonly_MooringDatabase_{MooringDatabase}.csv'
path1 = f'/Turbines_tol_{tolerance}_filter_designonly_MooringDatabase_{MooringDatabase}'
opt_baseline_file = path + '/BaselineOptimization/SNOPTlayoutnew_square_9T_iea.csv'

# Mooring Database inputs
Xdiff_file=path0+f'/MooringDatabase_{MooringDatabase}/WindSpeed_10'+ '/Xdiff_Surge.csv' # The positions achieved by each mooring design in each wind dir (x-axis)
Ydiff_file=path0+f'/MooringDatabase_{MooringDatabase}/WindSpeed_10'+ '/Xdiff_Sway.csv' # The positions achieved by each mooring design in each wind dir (y-axis)


shape='square'



Floris_perp= pd.read_csv(path+ '/SingleWIndDir_Opt'+opt_dyn_layout_file)
Floris_perp.columns=Floris_perp.columns.astype(str).astype(float)

layout_final=pd.read_csv(path+ '/FinalMooringWFDesign/' +outputfile,index_col=False)


fi = FlorisInterface('../Inputs/gch15MW.yaml')
fi.floris.flow_field.turbulence_intensity=0.06
cutin_ws=4
cutout_ws=25
# Windrose parameters
############################iea wind rose ############################

wd=[0, 22.5, 45, 67.5, 90, 112.5, 135, 157.5, 180, 202.5, 225, 247.5, 270, 292.5, 315, 337.5]
freq_d= [.025,  .024,  .029,  .036,.063,  .065,  .100,  .122,.063,  .038,  .039,  .083, .213,  .046,  .032,  .022]
k=np.zeros(16)+2.5
lam = np.zeros(16)+10

freq=np.transpose(np.array([freq_d]))
Windrose=make_wind_rose_from_weibull(wd, freq_d, wd=wd, ws=np.arange(4,26,1),
                                     weibull_k=k, weibull_lam=lam)


############################ wind rose plot ############################
legend_kwargs=dict(ncol=1, loc='best', fancybox=True,bbox_to_anchor=(0.4, 0.55, 0.45, 0.5))
mydpi=96
width=450.6/72*0.7
fig = plt.figure(figsize=(width,width),dpi=200)
ax = fig.add_subplot(projection='polar')
ws_right_edges = np.array([5, 10, 15, 20, 25])
color_array = sns.color_palette("rocket", n_colors = len(ws_right_edges))
ws_step = ws_right_edges[1] - ws_right_edges[0]
ws_labels = ["%d-%d m/s" % (w - ws_step, w) for w in ws_right_edges]

for wdir in wd:
    rects = []
    df_plot_sub = Windrose[Windrose.wd == wdir]
    for ws_idx, ws in enumerate(ws_right_edges[::-1]):
        plot_val = df_plot_sub[
            df_plot_sub.ws <= ws
        ].freq_val.sum()  # Get the sum of frequency up to this wind speed
        rects.append(
            ax.bar(
                np.radians(wdir),
                plot_val,
                width=1* np.radians(wd[1]-wd[0]),
                color=color_array[ws_idx],
                edgecolor="k",
            )
        )
ax.legend(reversed(rects), ws_labels, **legend_kwargs)
ax.set_theta_direction(-1)
ax.set_theta_offset(np.pi / 2.0)
ax.set_theta_zero_location("N")
ax.set_xticks(np.arange(0, 2*np.pi, np.pi/4))
ax.set_xticklabels(["N", "NE", "E", "SE", "S", "SW", "W", "NW"])
ax.set_rticks([0.04, 0.08, 0.12, 0.16])  
ax.set_rlabel_position(70)
# plt.savefig(r"/media/sf_Shared_folder"+'/PhD_Figs'+ '/Windrose_HornsRev.png',dpi=200,bbox_inches='tight')
plt.show()
########################################################


layoutdf=pd.read_csv(opt_baseline_file)

layout0=(np.array(layoutdf.x).tolist(),
        np.array(layoutdf.y).tolist())



wd_array = np.array(Windrose["wd"].unique(), dtype=float)
ws_array = np.array(Windrose["ws"].unique(), dtype=float)

wd_grid, ws_grid = np.meshgrid(wd_array, ws_array, indexing="ij")
freq_interp = NearestNDInterpolator(Windrose[["wd", "ws"]], Windrose["freq_val"])
freq = freq_interp(wd_grid, ws_grid)


if os.path.isfile(path0+'/AEP_nowake.csv') and os.path.isfile(path0+'/AEP_baseline.csv'):
    AEP_nowake = pd.read_csv(path0+'/AEP_nowake.csv')
    AEP_dir = pd.read_csv(path0+'/AEP_baseline.csv')

    
else:
    AEP_dir=pd.DataFrame()
    AEP_nowake=pd.DataFrame()
    for i in range(cutin_ws,26,1):
        TempWindrose=Windrose.loc[Windrose.ws==i,:].copy()
        new_wd=TempWindrose.wd.values
        new_freq=np.transpose(np.array([TempWindrose.freq_val.values]))
        new_ws=[i]
        for ww in range(len(new_wd)):  
            fi.reinitialize(layout_x=layout0[0],
                            layout_y=layout0[1],
                            wind_directions=[new_wd[ww]],
                            wind_speeds=new_ws) 
        
            AEP_dir.loc[i-cutin_ws,wd[ww]]=fi.get_farm_AEP([new_freq[ww]],freq_warn=0)
            AEP_nowake.loc[i-cutin_ws,wd[ww]]=fi.get_farm_AEP([new_freq[ww]],freq_warn=0,no_wake=True)
    
    AEP_dir.to_csv(path0+'/AEP_baseline.csv', index=False)
    AEP_nowake.to_csv(path0+'/AEP_nowake.csv', index=False)


AEP_initial = AEP_dir.sum().sum()
AEP_initial_nowake = AEP_nowake.sum().sum()
Wakeloss_baseline = pd.DataFrame((AEP_dir.values - AEP_nowake.values)/AEP_initial_nowake*100,columns=wd)


if os.path.isfile(path0+'/AEP_target.csv'):
    AEP_target = pd.read_csv(path0+'/AEP_target.csv')
    
else:
    AEP_target=pd.DataFrame()
    lx = pd.read_csv(path+ '/SingleWIndDir_Opt'+ '/OptimumDynamicfarmlayout_xaxis_square_9T_iea_TI_0.06.csv')
    ly = pd.read_csv(path+ '/SingleWIndDir_Opt'+ '/OptimumDynamicfarmlayout_yaxis_square_9T_iea_TI_0.06.csv')
    lx.columns = lx.columns.astype(float).astype(int)
    ly.columns = ly.columns.astype(float).astype(int)
    for i in range(cutin_ws,26,1):
        TempWindrose=Windrose.loc[Windrose.ws==i,:].copy()
        new_wd=TempWindrose.wd.values
        new_freq=np.transpose(np.array([TempWindrose.freq_val.values]))
        new_ws=[i]
        for ww in range(len(new_wd)):  
            w_dir= -new_wd[ww]+270
            w_dir%= 360
            fi.reinitialize(layout_x = lx.loc[:,int(myround(w_dir))],
                            layout_y = ly.loc[:,int(myround(w_dir))],
                            wind_directions=[new_wd[ww]],
                            wind_speeds=new_ws) 
        
            AEP_target.loc[i-cutin_ws,wd[ww]]=fi.get_farm_AEP([new_freq[ww]],freq_warn=0)
    
    AEP_target.to_csv(path0+'/AEP_target.csv', index=False)

AEP_total_target = AEP_target.sum().sum()
Wakeloss_target = pd.DataFrame((AEP_target.values - AEP_nowake.values)/AEP_initial_nowake*100,columns=wd)
Gaintarget = pd.DataFrame((AEP_target.values - AEP_dir.values)/AEP_dir.sum().sum()*100,columns=wd)

if os.path.isfile(path+"/FinalMooringWFDesign"+path1+'/AEP_final.csv'):
    AEP_new= pd.read_csv(path+"/FinalMooringWFDesign"+path1+'/AEP_final.csv')
else:
    Gain_Summary=pd.DataFrame(columns=['windspeed','gain'])
    AEP_new=pd.DataFrame(np.zeros((22,len(wd))),columns=wd)
    for i in range(cutin_ws,26,1):
        TempWindrose=Windrose.loc[Windrose.ws==i,:].copy()
        new_wd=TempWindrose.wd.values
        new_freq=np.transpose(np.array([TempWindrose.freq_val.values]))
        new_ws=[i]
        path2 =path+"/FinalMooringWFDesign"+path1+f"/WindSpeed_{i}" #if you want to compare exact directions
    
            
        Xdiff_file=path2+ '/Xdiff_Surge.csv' # The positions achieved by each mooring design in each wind dir (x-axis)
        Ydiff_file=path2+ '/Xdiff_Sway.csv' # The positions achieved by each mooring design in each wind dir (y-axis)
        Xdiff=pd.read_csv(Xdiff_file).round(decimals=2)
        Ydiff=pd.read_csv(Ydiff_file).round(decimals=2)
        
        Layoutdict=layoutfromMoorSystem(layout_final,wd,Xdiff,Ydiff)
    
        delta_E = delta_energy(new_wd,new_ws,new_freq,layout0,Layoutdict,0,AEP_initial,AEP_dir.loc[i-cutin_ws,:])  
    
        gain0=sum(delta_E['delta_AEP'])
        print('gain=')
        print(gain0)
        Gain_Summary.loc[i-cutin_ws,'windspeed']=new_ws[0]
        Gain_Summary.loc[i-cutin_ws,'gain']=gain0
        AEP_new.loc[i-cutin_ws,:] = delta_E['AEP_new'].values
        
    Total_gain=sum(Gain_Summary['gain'])
    AEP_new.to_csv(path+"/FinalMooringWFDesign"+path1+'/AEP_final.csv', index=False)

Wakeloss_new = pd.DataFrame((AEP_new.values - AEP_nowake.values)/AEP_initial_nowake*100,columns=wd)
Gainmatrix = pd.DataFrame((AEP_new.values - AEP_dir.values)/AEP_dir.sum().sum()*100,columns=wd)

############## new part ##############
if os.path.isfile(path+"/FinalMooringWFDesign"+path1+'/AEP_final_withmotions.csv'):
    AEP_new_updated= pd.read_csv(path+"/FinalMooringWFDesign"+path1+'/AEP_final_withmotions.csv')
else:
    Gain_Summary=pd.DataFrame(columns=['windspeed','gain'])
    AEP_new_updated=pd.DataFrame(np.zeros((22,len(wd))),columns=wd)
    for i in range(cutin_ws,26,1):
        TempWindrose=Windrose.loc[Windrose.ws==i,:].copy()
        new_wd=TempWindrose.wd.values
        new_freq=np.transpose(np.array([TempWindrose.freq_val.values]))
        new_ws=[i]
        path2 =path+"/FinalMooringWFDesign"+path1+f"/WindSpeed_{i}" #if you want to compare exact directions
    
            
        Xdiff_file=path2+ '/Xdiff_Surge.csv' # The positions achieved by each mooring design in each wind dir (x-axis)
        Ydiff_file=path2+ '/Xdiff_Sway.csv' # The positions achieved by each mooring design in each wind dir (y-axis)
        Sway_file=path2+ '/X_Sway.csv' # The positions achieved by each mooring design in each wind dir (y-axis)
        Xdiff=pd.read_csv(Xdiff_file).round(decimals=2)
        Ydiff=pd.read_csv(Ydiff_file).round(decimals=2)
        Sway=pd.read_csv(Sway_file).round(decimals=2)
        Xdiff_new = Xdiff.copy()
        Ydiff_new = Ydiff.copy()
        Sway_new = Sway.copy()
        Layoutdict = layoutfromMoorSystem(layout_final,wd,Xdiff,Ydiff)    
        T_velocity = round(turbine_velocity(new_wd,new_ws,new_freq,layout0,Layoutdict),4)
        vnew = T_velocity * 0
        vprev = T_velocity * 0 + 10
        counter = 0
        constraint =10
        while constraint !=0 and i < 15:
            print(f'number of iterations: {counter}')
            vprev = T_velocity
            for ww in range(len(new_wd)):  
                w_dir= -new_wd[ww]+270
                w_dir%= 360
    
                for t in range(len(Xdiff)):
                    ws_low= int(np.floor(T_velocity.loc[t,w_dir]))
                    ws_high = int(np.ceil(T_velocity.loc[t,w_dir]))
                    if ws_low<4 or ws_high <4:
                        ws_low = ws_high = 4
                        Xdiff_temp = pd.read_csv( path+"/FinalMooringWFDesign"+path1+f"/WindSpeed_{ws_low}" + '/Xdiff_Surge.csv').round(decimals=2)
                        Ydiff_temp = pd.read_csv( path+"/FinalMooringWFDesign"+path1+f"/WindSpeed_{ws_low}" + '/Xdiff_Sway.csv' ).round(decimals=2)
                        Sway_temp = pd.read_csv( path+"/FinalMooringWFDesign"+path1+f"/WindSpeed_{ws_low}" + '/X_Sway.csv' ).round(decimals=2)
                        w_d= np.round(w_dir/5)*5
                        Xdiff_new.loc[t,str(int(w_d))] = Xdiff_temp.loc[t,str(int(w_d))] 
                        Ydiff_new.loc[t,str(int(w_d))] = Ydiff_temp.loc[t,str(int(w_d))] 
                        Sway_new.loc[t,str(int(w_d))] = Sway_temp.loc[t,str(int(w_d))] 
    
                    elif ws_low==ws_high:
                        Xdiff_temp = pd.read_csv( path+"/FinalMooringWFDesign"+path1+f"/WindSpeed_{ws_low}" + '/Xdiff_Surge.csv').round(decimals=2)
                        Ydiff_temp = pd.read_csv( path+"/FinalMooringWFDesign"+path1+f"/WindSpeed_{ws_low}" + '/Xdiff_Sway.csv' ).round(decimals=2)
                        Sway_temp = pd.read_csv( path+"/FinalMooringWFDesign"+path1+f"/WindSpeed_{ws_low}" + '/X_Sway.csv' ).round(decimals=2)
                        w_d= np.round(w_dir/5)*5
                        Xdiff_new.loc[t,str(int(w_d))] = Xdiff_temp.loc[t,str(int(w_d))] 
                        Ydiff_new.loc[t,str(int(w_d))] = Ydiff_temp.loc[t,str(int(w_d))] 
                        Sway_new.loc[t,str(int(w_d))] = Sway_temp.loc[t,str(int(w_d))] 
                    else:
                        Xdiff_low = pd.read_csv( path+"/FinalMooringWFDesign"+path1+f"/WindSpeed_{ws_low}" + '/Xdiff_Surge.csv').round(decimals=2)
                        Ydiff_low = pd.read_csv( path+"/FinalMooringWFDesign"+path1+f"/WindSpeed_{ws_low}" + '/Xdiff_Sway.csv' ).round(decimals=2)
                        Sway_low = pd.read_csv( path+"/FinalMooringWFDesign"+path1+f"/WindSpeed_{ws_low}" + '/X_Sway.csv' ).round(decimals=2)
                        Xdiff_high = pd.read_csv( path+"/FinalMooringWFDesign"+path1+f"/WindSpeed_{ws_high}" + '/Xdiff_Surge.csv').round(decimals=2)
                        Ydiff_high = pd.read_csv( path+"/FinalMooringWFDesign"+path1+f"/WindSpeed_{ws_high}" + '/Xdiff_Sway.csv' ).round(decimals=2)
                        Sway_high = pd.read_csv( path+"/FinalMooringWFDesign"+path1+f"/WindSpeed_{ws_high}" + '/X_Sway.csv' ).round(decimals=2)
                        w_d= np.round(w_dir/5)*5
                        Xdiff_new.loc[t,str(int(w_d))] = Xdiff_low.loc[t,str(int(w_d))] + ((T_velocity.loc[t,w_dir]-ws_low)*(Xdiff_high.loc[t,str(int(w_d))]-Xdiff_low.loc[t,str(int(w_d))])/(ws_high-ws_low))
                        Ydiff_new.loc[t,str(int(w_d))] = Ydiff_low.loc[t,str(int(w_d))] + ((T_velocity.loc[t,w_dir]-ws_low)*(Ydiff_high.loc[t,str(int(w_d))]-Ydiff_low.loc[t,str(int(w_d))])/(ws_high-ws_low))
                        Sway_new.loc[t,str(int(w_d))] = Sway_low.loc[t,str(int(w_d))] + ((T_velocity.loc[t,w_dir]-ws_low)*(Sway_high.loc[t,str(int(w_d))]-Sway_low.loc[t,str(int(w_d))])/(ws_high-ws_low))
            Layoutdict = layoutfromMoorSystem(layout_final,wd,Xdiff_new,Ydiff_new)    
            T_velocity = round(turbine_velocity(new_wd,new_ws,new_freq,layout0,Layoutdict),4)
            vnew = T_velocity
            counter = counter + 1 
            constraint = round(sum(sum(abs(vnew.values-vprev.values))),1)
    
        
        if i <= 14:
            Layoutdict = layoutfromMoorSystem(layout_final,wd,Xdiff_new,Ydiff_new)   
        delta_E = delta_energy(new_wd,new_ws,new_freq,layout0,Layoutdict,0,AEP_initial,AEP_dir.loc[i-cutin_ws,:]) 
        T_velocity = round(turbine_velocity(new_wd,new_ws,new_freq,layout0,Layoutdict),4)
        if i ==10:
            T_velocity.to_csv(path+"/FinalMooringWFDesign"+path1+f'/floris_rotorspeed_atws_{i}.csv', index=False)
    
        gain0=sum(delta_E['delta_AEP'])
        print('gain=')
        print(gain0)
        Gain_Summary.loc[i-cutin_ws,'windspeed']=new_ws[0]
        Gain_Summary.loc[i-cutin_ws,'gain']=gain0
        AEP_new_updated.loc[i-cutin_ws,:] = delta_E['AEP_new'].values
        
    Total_gain=sum(Gain_Summary['gain'])
    AEP_new_updated.to_csv(path+"/FinalMooringWFDesign"+path1+'/AEP_final_withmotions.csv', index=False)

Wakeloss_new_updated = pd.DataFrame((AEP_new_updated.values - AEP_nowake.values)/AEP_initial_nowake*100,columns=wd)
Gainmatrix_updated = pd.DataFrame((AEP_new_updated.values - AEP_dir.values)/AEP_dir.sum().sum()*100,columns=wd)

