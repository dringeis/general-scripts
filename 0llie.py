#! /usr/bin/env python

from math import ceil

### Parameters

## Simulation Domain sizes
X=100 #Domain size in x
Y=250 #Domain size in y
Lx=[2,100] # limits in the search of cases for Nx
# Ly=[20,1500] # limits in the search of cases for Ny
Ly=Lx # limits in the search of cases for Ny

## Parameters for the parallelism
Nmax=128 # maximum number of node you want to consider
Nmin=1 # minimum number of node you want to consider

Pnmin=24 # Minimum number of job per node

even_only=False # Only show the cases with even distribution of tasks

### Computations :

print 'Domain : ',X,'x',Y
print ' '
print 'Solutions to optimize the use of ollie :'
print '   Nx and Ny, the dimension of MPI divisions in X and Y'
print '   Px and Py, the numbers of cores in X and Y'
print '   Ptot and N, the total number of core used and the total number of nodes'
print '   I , the number of unused cores in this solution'
print ' '

Np=36 # number of cores per nodes (ollie)

f=1 # factor to expand the search with for unused cores superior to one node (even runs only)
for Ip in xrange(f*Np+1):
  for Nx in xrange(Lx[0],Lx[1]+1):
    if X%Nx==0:
      Px=X/Nx
      for Ny in xrange(Ly[0],Ly[1]+1):
        if Y%Ny==0:
          Py=Y/Ny
          Ptot=float(Px*Py)
          N=int(ceil(Ptot/Np))
          for Ni in xrange(max(N,Nmin),N+f+1):
            I=Ni*Np-int(Ptot)
            Pn=int(Ptot/Ni)
            In=int(I/Ni)
            if I==Ip and Ni >= Nmin and Ni <= Nmax and Pn >= Pnmin:
              if I%Ni==0 :
                  print 'Nx :', Nx, ',Ny :',Ny, ',Px :',Px, ',Py :',Py, ',Ptot :',Ptot, ',N :', Ni ,',I :',I ,',RQ :',I/Ptot,' >> good solution :' , Pn , ' running jobs per node ,', In , 'idle cores per node'
              elif I<Np :
                if even_only!=True :
                  print 'Nx :', Nx, ',Ny :',Ny, ',Px :',Px, ',Py :',Py, ',Ptot :',Ptot, ',N :', Ni, ',I :',I ,',RQ :',I/Ptot
