# general-scripts
Various scripts I write to work on MITgcm sea ice simulations.
Without waranty 

## List of files and descriptions :
  - *plot_data.py*: Various functions to open, read, select variables in order to plot seaice variable base on the xmitgcm python package. It also includes routines to plot the yield curve.
  - *ploti.py*: Interactive ploting code, using the routines from *plot_data.py*. Read the data and data.diagnostics to get the timetstep the diagnostics files to open.
  - *yield_curve.py*: Plotting the yield curve from a generic set of strain rates, convienient to test your rheology changes before implementing them in MITgcm
  - *0llie.py*: This simple piece of code is made to explore the different mpi possibilities depending in the size of your domain, in order to maximize the usage of the super computer. It is designed for AWI super computer, but can be easily tuned for other supercomputers.
  
