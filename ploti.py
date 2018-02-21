#! /usr/bin/env python

# import std librairies
import sys
import os
import getopt

import subprocess

import warnings
warnings.simplefilter("ignore")

# loading plot_data routines
from plot_data import data_dt,data_diags,xm_open,xm_choose,xm_selec,xm_ploting

idiag = ''
var = ''
timestep = ''

try:
    opts, args = getopt.getopt(sys.argv[1:],"hi:v:t:",["help","ifile=","var=","tstp="])
except getopt.GetoptError:
    print 'Direct use : ploti.py -i <inputfile> -v <variable> -t <timestep>'
for opt, arg in opts:
    if opt in ('-h','--help'):
       print 'ploti.py -i <inputfile> -v <variable> -t <timestep>'
       sys.exit()
    elif opt in ("-i", "--ifile"):
        idiag = arg
    elif opt in ("-v", "--var"):
        var = arg
    elif opt in ("-t", "--tstp"):
        timestep = arg

if len(idiag)==0:
    diagfiles=data_diags()
    if len(diagfiles)==0:
        os.system('ls *.*.data')
        idiag=str(raw_input("Diagnostic file name to open :"))
    else :
        print "Opening DiagFiles : ",diagfiles
        idiag=diagfiles

dt=data_dt()
print "extracting dt from data : dt = ",dt

data=xm_open(idiag,dt)
print "Size of the Domain X x Y x Z , dt : ",data.dims['XC'],' x ',data.dims['YC'],' x ',data.dims['Z'],' , ',dt
print "Number of timesteps available :",data.dims['time']-1
print data.data_vars

if len(var)==0 or len(timestep)==0:
    var_name, tstp=xm_choose(data,dt,var,timestep)
else:
    var_name=var
    tstp=int(timestep)

#var_t,ticks,var_max,var_min,varname,var_time,var_tstp,loga=xm_selec(data,var_name,tstp,dt)
#xm_ploting(var_t,ticks,var_max,var_min,varname,var_time,var_tstp,loga=loga)

while var_name != 'None':

    if var_name != 'None' :

        while tstp >= 0 :

            if tstp >= 0 or var_name != 'None' :

                var_t,ticks,var_max,var_min,varname,var_time,var_tstp,loga=xm_selec(data,var_name,tstp,dt)
                xm_ploting(var_t,ticks,var_max,var_min,varname,var_time,var_tstp,loga=loga)

            var_name, tstp=xm_choose(data,dt,var_name,[])

    var_name, tstp=xm_choose(data,dt,[],[])




