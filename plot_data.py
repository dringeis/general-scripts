#! /usr/bin/env python

# import std librairies
import sys
import os
import matplotlib
#matplotlib.use("Agg")
matplotlib.use("TkAgg")
import numpy as np
import pylab as plt
from pprint import pprint
import shutil
from copy import copy,deepcopy

from matplotlib.colors import SymLogNorm
from matplotlib import ticker



# Using LaTeX in figures
plt.rc('text', usetex=True)
plt.rc('font', family='sans')

# adding my directoty with python files
sys.path.append(os.path.abspath("~/init_code/"))

# import librairy to read output file from MITgcm
from python_configuration import *

# import the awesome new colormaps
from colormaps import viridis


'''
#### xmitgcm package

## to install:

pip install xarray --user
pip install dask --user
pip install xmitgcm --user


'''

from xmitgcm import open_mdsdataset


def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'$10^{{{}}}$'.format(b)


def data_dt():
    searchfile = open("data", "r")
    for line in searchfile:
        if "deltaT" in line: deltaTs=line
    searchfile.close()

    ll=deltaTs.split("=")
    dt=ll[-1].replace(',','')

    return float(dt)



def data_diags():
    searchfile = open("data.diagnostics", "r")
    files=[]
    for line in searchfile:
        if "fileName" in line and '#' not in line:
            name=line.split("=")[-1].replace("'",'').replace(',','').replace(' ','').replace('\n','')
            files.append(name)
    searchfile.close()

    return files



def xm_open(name_file,dt):
    data_dir = './'
    data = open_mdsdataset(data_dir,prefix=name_file,delta_t=dt)

    return data



def xm_choose(data,dt,var,tstp):

    if len(var)==0:
        print('Enter None instead of variable name to go out of plotting')
        var=str(raw_input("Desired Variable: "))

    if len(tstp)==0:
        print('Enter negative timestep to go out of ploting current Variable')
        tstp=int(raw_input("Desired Timestep :"))


    return var, tstp



def xm_selec(data,varname,tstep,dt=0.0):

    var_t=data[varname][tstep]

    ticks=None

    if varname=='SIpress':
        loga=False
        Heff_t=data.SIheff[tstep]
        press_real=var_t/Heff_t
        var_t=press_real

    if varname=='SIheff':
        loga=False
        var_max=1.5
        var_min=0.0
    elif varname=='SIarea':
        loga=False
        var_max=1.0
        var_min=0.0
    elif varname=='SIshear':
        loga=True
        var_max=0.01
        var_min=0.0
    elif varname=='SIdiv':
        loga=True
        var_max=0.020
        var_min=-0.020
    elif varname=='SIstress':
        loga=False
        var_max=0.1
        var_min=0.0
    else :
        loga=False
        var_max=var_t.values.max()
        var_min=var_t.values.min()

    n=6
    if var_min==0.0: var_min=1e-7
    ticks=np.logspace(np.log10(var_min),np.log10(var_max),n)

    var_tstp=data.coords['iter'][tstep].values
    var_time=data.coords['time'][tstep].values

    return var_t,ticks,var_max,var_min,varname,var_time,var_tstp,loga



def xm_ploting(var_t,ticks,var_max,var_min,varname,var_time,var_tstp,save=False,savefile=False,loga=False,title=True):

    plt.figure()

    if loga==False:
        fig=var_t.plot(vmin=var_min,vmax=var_max)
    elif loga==True:
        fig=var_t.plot(norm=SymLogNorm(vmin=var_min,vmax=var_max,linthresh=1e-9,linscale=1e-9),cbar_kwargs={'ticks':ticks,'format':ticker.FuncFormatter(fmt)})


    if title==True:
        if var_time==0.0 and var_tstp!=0 :
            plt.title(r"\textbf{Variable} "+str(varname)+r" \textbf{at }"+str(var_tstp).zfill(9)   +' timesteps')
        else :
            plt.title(r"\textbf{Variable} "+str(varname)+r" \textbf{at }"+str(var_time).zfill(6)   +' s')
    else:
        plt.title(' ')

    plt.axes().set_aspect('equal', 'datalim')

    plt.xlabel(r'x [m]')
    plt.ylabel(r'y [m]')

    if save == True :
        if savefile==False:
            plt.savefig(str(str(varname)+"_"+str(var_tstp).zfill(9)+".png"),bbox_inches='tight',dpi=300) # saving figure
        else :
            plt.savefig(savefile,bbox_inches='tight',dpi=300) # saving figure
    elif save == False :
        plt.show()

    return None




'''
Function to plot data from MITgcm using the diagnostics output and the xmitgcm package
xmitgcm package : https://xmitgcm.readthedocs.io/en/latest/

ONLY FOR SEA ICE DATA : no depth support, 3D (T,Y,X) only, not 4D (T,Z,Y,X)

input :
    - name_diag - string  - name of the diagnostics file - ex : 'iceDiag'
    - varname - string- name of the variable to plot and/or save  - ex : 'SIheff'
    - tstep - int - the timestep to plot

    - options :
        - loga - boolean- if you want the variable to be ploted with log colorbar
            def : False
        - save - boolean - if you want to save the picture and not show it.
            def : False
        - dt - double - timestep duration (in seconds)
            def : 0.0 (print timesteps)

output :
    - None but plot the variable in directories with the shape :
        'plots_VARIABLE_TIMESTEP-BEGIN_TIMESTEP-END/'
    - with name of the shape :
        'VARNAME_TIMESTEPNUMBER.png'

'''
def plot_xm(name_file,varname,tstep,loga=False,save=False,dt=0.0,title=True,savefile=False):

    data_dir = './'
    data = open_mdsdataset(data_dir,prefix=[name_file])

    var_t=data[varname][tstep]

    ticks=None

    if varname=='SIpress':
        Heff_t=data.SIheff[tstep]
        press_real=var_t/Heff_t
        var_t=press_real


    if varname=='SIheff':
        var_max=1.5
        var_min=0.0
    elif varname=='SIarea':
        var_max=1.0
        var_min=0.0
    elif varname=='SIshear':
        # ticks=[0.00,0.001,0.010,0.020,0.050,0.16]
        var_max=0.01
        var_min=0.0
    elif varname=='SIdiv':
        var_max=0.020
        var_min=-0.020
    elif varname=='SIstress':
        var_max=0.1
        var_min=0.0
    else :
        var_max=var_t.values.max()
        var_min=var_t.values.min()

    n=6
    if var_min==0.0: var_min=1e-7
    ticks=np.logspace(np.log10(var_min),np.log10(var_max),n)

    var_tstp=data.coords['time'][tstep].values
    var_time=var_tstp*dt

    plot_xm_plotting(var_t,ticks,var_max,var_min,varname,var_time,var_tstp,save=save,savefile=savefile)

    ## replaced by a function above
    # plt.figure()

    # if loga==False:
    #     fig=var_t.plot(vmin=var_min,vmax=var_max)
    # elif loga==True:
    #     fig=var_t.plot(norm=SymLogNorm(vmin=var_min,vmax=var_max,linthresh=1e-9,linscale=1e-9),cbar_kwargs={'ticks':ticks,'format':ticker.FuncFormatter(fmt)})


    # if title==True:
    #     if dt==0.0 :
    #         plt.title(r"\textbf{Variable} "+str(varname)+r" \textbf{at }"+str(var_tstp).zfill(9)   +' timesteps')
    #     else :
    #         plt.title(r"\textbf{Variable} "+str(varname)+r" \textbf{at }"+str(var_time).zfill(6)   +' s')
    # else:
    #     plt.title(' ')

    # plt.axes().set_aspect('equal', 'datalim')

    # plt.xlabel(r'x [m]')
    # plt.ylabel(r'y [m]')

    # if save == True :
    #     if savefile==False:
    #         plt.savefig(str(str(varname)+"_"+str(var_stsp).zfill(9)+".png"),bbox_inches='tight',dpi=300) # saving figure
    #     else :
    #         plt.savefig(savefile,bbox_inches='tight',dpi=300) # saving figure
    # elif save == False :
    #     plt.show()

    return None





'''
Function to plot data from MITgcm using the diagnostics output and the xmitgcm package
xmitgcm package : https://xmitgcm.readthedocs.io/en/latest/

ONLY FOR SEA ICE DATA : no depth support, 3D (T,Y,X) only, not 4D (T,Z,Y,X)

This function is called by plot_var_xm below

input :
    - data - xmitgcm object (dasks array in xarray structure)  - xmitgcm output - data
    - varname - string - names of the variable to plot and save  - ex : 'SIheff'
    - tstep - int - tstep to plot in [0,n-1] - ex : 10
    - loga - boolean - if you want the variable to be ploted with log colorbar - ex : False
    - savedir - string - saving directory
    - var_min - double - Lower limit for colorbar
    - var_max - double - Upper limit for colorbar
    - dt - double - duration of timestep

output :
    - None but plot the variable in directories with the shape :
        'plots_VARIABLE_TIMETSTEPBEGIN_TIMESTEPEND/'
    - with name of the shape :
        'VARNAME_TIMESTEPNUMBER.png'
'''

def plot_xm_forloop(data,varname,tstep,loga,savedir,var_min,var_max,dpi,dt):

    var=data[varname]
    var_t=data[varname][tstep]

    if varname=='SIpress':
        Heff_t=data['SIheff'][tstep]
        press_real=var_t/Heff_t
        var_t=press_real

    plt.figure()
    if loga==False :
        var_t.plot(vmin=var_min,vmax=var_max)
    if loga==True:
        var_t.plot(vmin=var_min,vmax=var_max,norm=SymLogNorm(vmin=var_min,vmax=var_max,linthresh=1e-8,linscale=1e-8))

    if dt==0.0:
        plt.title(r"\textbf{Variable} "+str(varname)+r" \textbf{at }"+str(data.coords['time'][tstep].values).zfill(9)   +' timesteps')
    else:
        plt.title(r"\textbf{Variable} "+str(varname)+r" \textbf{at }"+str(data.coords['time'][tstep].values*dt).zfill(9)   +' s')

    plt.axes().set_aspect('equal', 'datalim')

    plt.xlabel('x [m]')
    plt.ylabel('y [m]')

    plt.savefig(savedir+str(str(varname)+"_"+str(tstep).zfill(9)+".png"),bbox_inches='tight',dpi=dpi)
    plt.close('all')

    return None

'''
Function to plot all a group of data from MITgcm using the diagnostics output and the xmitgcm package
xmitgcm package : https://xmitgcm.readthedocs.io/en/latest/

ONLY FOR SEA ICE DATA : no depth support, 3D (T,Y,X) only, not 4D (T,Z,Y,X)

Is calling the function 'plot_xm_forloop' above.

input :
    - name_file - string  - name of the diagnostics file - ex : 'iceDiag'
    - varnames - string-array - names of the variable to plot and save  - ex : ['SIheff','SIpress']

    - Options :
        - loga - boolean-array - if you want the variable to be ploted   - ex : [False, True]
            def : all false
        - beg - int - first timestep to be ploted
            def : 0
        - end - int - last time step to be ploted
            def : -1 (last)
        - step - int - step to plot the figures
            def : 1
        - dpi - int - Quality of output images
            def : None (matplotlib default)
        - dt - double - timestep duration (in seconds)
            def : 0.0 (print timesteps)

output :
    - None but plot the variable in directories with the shape :
        'plots_VARIABLE_TIMETSTEPBEGIN_TIMESTEPEND/'
    - with name of the shape :
        'VARNAME_TIMESTEPNUMBER.png'

'''
def plot_var_xm(name_file,varnames,loga=0.,beg=0,end=-1,step=1,dpi=None,dt=0.0):

    data_dir = './'
    data = open_mdsdataset(data_dir,prefix=[name_file])

    if loga==0. :
        loga=np.full((len(varnames)), False, dtype=bool)

    if beg==0 and end ==-1 and step==1 :
        nbr=data.dims['time']
        times = xrange(nbr)

    else :
        var_cho=np.array(var[beg:end+1:step,:])

    for ivar in range(len(varnames)):
        name_var=varnames[ivar]

        var_dir=str('./plots_'+name_var+'_'+str(data.coords['iter'][beg].values)+'_'+str(data.coords['iter'][end].values))
        var_dirp=str(var_dir+'/')

        if os.path.exists(var_dir):
            shutil.rmtree(var_dir)
        os.mkdir(var_dir)

        plotb=False

        var=data[name_var]
        if name_var=='SIheff':
            var_max=2.0
            var_min=0.0
            plotb=True
        elif name_var=='SIarea':
            var_max=1.0
            var_min=0.0
            plotb=True
        elif name_var=='SIshear':
            var_max=0.020
            var_min=0.0
            plotb=True
        elif name_var=='SIdiv':
            var_max=0.020
            var_min=-0.020
            plotb=True
        elif name_var=='SIstress':
            var_max=0.1
            var_min=0.0
            plotb=True
        else:
            print 'Enter vmin and vmax values for this variable!'
            plotb=False

        if plotb :
            for tstep in times:
                plot_xm_forloop(data,name_var,tstep,loga[ivar],var_dirp,var_min,var_max,dpi,dt)

    return None



from matplotlib.patches import Ellipse

'''
Function to plot sea ice VP yield curve from MITgcm simulations using the diagnostics output and the xmitgcm package
xmitgcm package : https://xmitgcm.readthedocs.io/en/latest/

ONLY FOR SEA ICE DATA : no depth support, 3D (T,Y,X) only, not 4D (T,Z,Y,X)

input :
    - name_diag - string  - name of the diagnostics file - ex : 'iceDiag'
    - varname - string-array - names of the variable to plot and save  - ex : 'SIheff'
    - tstep - int - the timestep to plot

    - options :
        - loga - boolean- if you want the variable to be ploted with log colorbar
            def : False
        - save - boolean - if you want to save the picture and not show it.
            def : False


output :
    - None but plot the variable in directories with the shape :
        'plots_VARIABLE_TIMESTEP-BEGIN_TIMESTEP-END/'
    - with name of the shape :
        'VARNAME_TIMESTEPNUMBER.png'

'''
def plot_yield(name_file,data_dir = './',save=False,e=2,f=1.0,t=0.0,tstp=-1):

    data = open_mdsdataset(data_dir,prefix=[name_file],chunks=None)

    sig1=data['SIsig1'][tstp]
    sig2=data['SIsig2'][tstp]

    sig1f=sig1.values.flatten()
    sig2f=sig2.values.flatten()

    print len(sig1f)

    plt.figure(1)
    ax = plt.gca()

    elli=Ellipse(xy=((-f+t)/2,(-f+t)/2), width=(f+t)/e*np.sqrt(2), height=(f+t)*np.sqrt(2), angle=-45,edgecolor='b', fc='None', lw=0.5)

    ax.add_patch(elli)
    plt.plot(sig1f,sig2f,'r.',markersize=1)

    xlim=ax.get_xlim()
    ylim=ax.get_ylim()

    ax.set(xlim=xlim, ylim=ylim)

    cordx=[min(xlim[0],ylim[0]),min(xlim[1],ylim[1])]
    cordy=[min(xlim[0],ylim[0]),min(xlim[1],ylim[1])]
    plt.plot(cordx,cordy,'-k',lw=1)

    coradx=[-min(-xlim[0],ylim[1]),min(xlim[1],-ylim[0])]
    corady=[min(-xlim[0],ylim[1]),-min(xlim[1],-ylim[0])]
    plt.plot(coradx,corady,'-k',lw=1)

    plt.axvline(linewidth=1, color='k')
    plt.axhline(linewidth=1, color='k')

    plt.title( r" Yield curve ")
    plt.axes().set_aspect('equal')
    plt.grid()

    if save == True :
        plt.savefig("yield_curve.png",bbox_inches='tight',dpi=200)

    plt.show()

#    plt.close('all')

    return None


import matplotlib.animation as manimation
matplotlib.use("Agg")

'''
Function to plot sea ice VP yield curve animation from MITgcm simulations using the diagnostics output and the xmitgcm package
xmitgcm package : https://xmitgcm.readthedocs.io/en/latest/

ONLY FOR SEA ICE DATA : no depth support, 3D (T,Y,X) only, not 4D (T,Z,Y,X)

input :
    - name_diag - string  - name of the diagnostics file - ex : 'iceDiag'
    - varname - string-array - names of the variable to plot and save  - ex : 'SIheff'
    - tstep - int - the timestep to plot

    - options :
        - loga - boolean- if you want the variable to be ploted with log colorbar
            def : False
        - save - boolean - if you want to save the picture and not show it.
            def : False


output :
    - None but plot the variable in directories with the shape :
        'plots_VARIABLE_TIMESTEP-BEGIN_TIMESTEP-END/'
    - with name of the shape :
        'VARNAME_TIMESTEPNUMBER.png'

'''
def plot_yield_vid(name_file,data_dir = './',e=2,f=1.0,t=0.0):

    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title='yield', artist='dringeis')
    writer = FFMpegWriter(fps=30, metadata=metadata,bitrate=-1)

    data = open_mdsdataset(data_dir,prefix=[name_file],chunks=None)

    sig1=data['SIsig1']
    sig2=data['SIsig2']

    fig=plt.figure()
    ax = fig.gca()

    elli=Ellipse(xy=((-f+t)/2,(-f+t)/2), width=(f+t)/e*np.sqrt(2), height=(f+t)*np.sqrt(2), angle=-45,edgecolor='b', fc='None', lw=0.5)

    ax.add_patch(elli)
    im, = plt.plot([],[],'r.',markersize=1)

    xlim=ax.get_xlim()
    ylim=ax.get_ylim()

    ax.set(xlim=xlim, ylim=ylim)

    cordx=[min(xlim[0],ylim[0]),min(xlim[1],ylim[1])]
    cordy=[min(xlim[0],ylim[0]),min(xlim[1],ylim[1])]
    plt.plot(cordx,cordy,'-k',lw=0.5)

    coradx=[-min(-xlim[0],ylim[1]),min(xlim[1],-ylim[0])]
    corady=[min(-xlim[0],ylim[1]),-min(xlim[1],-ylim[0])]
    plt.plot(coradx,corady,'-k',lw=0.5)

    plt.axvline(linewidth=0.5, color='k')
    plt.axhline(linewidth=0.5, color='k')

    plt.title( r" Yield curve ")
    plt.axes().set_aspect('equal')
    plt.grid()

    with writer.saving(fig, "yield.mp4", 200):
        for t in xrange(data.dims['time']):
            sig1f=sig1[t,:,:].values.flatten()
            sig2f=sig2[t,:,:].values.flatten()
            im.set_data(sig1f, sig2f)
            writer.grab_frame()

    return None

'''
Function to plot comparison data from 2 MITgcm experiment using the diagnostics output and the xmitgcm package
xmitgcm package : https://xmitgcm.readthedocs.io/en/latest/

ONLY FOR SEA ICE DATA : no depth support, 3D (T,Y,X) only, not 4D (T,Z,Y,X)

This function is called by plot_comp below

input :
    - data - xmitgcm object (dasks array in xarray structure)  - xmitgcm output - data
    - varname - string - names of the variable to plot and save  - ex : 'SIheff'
    - tstep - int - tstep to plot in [0,n-1] - ex : 10
    - loga - boolean - if you want the variable to be ploted with log colorbar - ex : False
    - savedir - string - saving directory

output :
    - None but plot the variable in directories with the shape :
        'plots_comp_VARIABLE_TIMETSTEPBEGIN_TIMESTEPEND/'
    - with name of the shape :
        'VARNAME_TIMESTEPNUMBER.png'
'''
def plot_comp_forloop(data1,data2,varname,tstep,loga,savedir,var_min,var_max,dpi):

    data_var_t1=data1[varname][tstep]
    data_var_t2=data2[varname][tstep]
    var_t=data_var_t2/data_var_t1

    plt.figure()
    if loga==False :
        var_t.plot(vmin=var_min,vmax=var_max)
    if loga==True:
        var_t.plot(vmin=var_min,vmax=var_max,norm=SymLogNorm(vmin=var_min,vmax=var_max,linthresh=1e-8,linscale=1e-8))

    plt.title(r"\textbf{Variable} "+str(varname)+r" \textbf{at }"+str(data1.coords['time'][tstep].values).zfill(9)   +' s')
    plt.axes().set_aspect('equal', 'datalim')
    plt.savefig(savedir+str(str(varname)+"_"+str(tstep).zfill(9)+".png"),bbox_inches='tight',dpi=dpi)
    plt.close('all')

    return None


'''
Function to plot all a group of comparion data from 2 MITgcm  experiment using the diagnostics output and the xmitgcm package
xmitgcm package : https://xmitgcm.readthedocs.io/en/latest/

ONLY FOR SEA ICE DATA : no depth support, 3D (T,Y,X) only, not 4D (T,Z,Y,X)

Is calling the function 'plot_comp_forloop' above.

input :
    - name_file - string  - name of the diagnostics file - ex : 'iceDiag'
    - dirs - string-array - names of the ditrectories to compare, plot and save  - ex : ['run00','run01']
       !! Is doing dir1/dir2 !!
    - varnames - string-array - names of the variable to plot and save  - ex : ['SIheff','SIpress']

    - Options :
        - loga - boolean-array - if you want the variable to be ploted   - ex : [False, True]
            def : all false
        - beg - int - first timestep to be ploted
            def : 0
        - end - int - last time step to be ploted
            def : -1 (last)
        - step - int - step to plot the figures
            def : 1

output :
    - None but plot the variable in directories with the shape :
        'plots_VARIABLE_TIMETSTEPBEGIN_TIMESTEPEND/'
    - with name of the shape :
        'VARNAME_TIMESTEPNUMBER.png'


'''
def plot_comp(name_file,dirs,varnames,loga=0.,beg=0,end=-1,step=1,dpi=100):

    data_dir1 = dirs[0]
    data1 = open_mdsdataset(data_dir1,prefix=[name_file])

    data_dir2 = dirs[1]
    data2 = open_mdsdataset(data_dir2,prefix=[name_file])

    if loga==0. :
        loga=np.full((len(varnames)), False, dtype=bool)

    if beg==0 and end ==-1 and step==1 :
        nbr=data1.dims['time']
        times = xrange(nbr)

    else :
        var_cho=np.array(var[beg:end+1:step,:])

    for ivar in range(len(varnames)):
        name_var=varnames[ivar]

        var_dir=str('./plots_comp_'+name_var+'_'+str(data1.coords['iter'][beg].values)+'_'+str(data1.coords['iter'][end].values))
        var_dirp=str(var_dir+'/')

        if os.path.exists(var_dir):
            shutil.rmtree(var_dir)
        os.mkdir(var_dir)

        var_max=0.9
        var_min=1.1

        for tstep in times:
            plot_comp_forloop(data1,data2,name_var,tstep,loga[ivar],var_dirp,var_min,var_max,dpi)

    return None

