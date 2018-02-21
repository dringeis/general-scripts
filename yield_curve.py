#!/usr/bin python
import numpy as np
from smoothabs import smoothmin_good
from matplotlib.patches import Ellipse
import pylab as plt

### Creating fake data
eg=np.mgrid[-100:100:200j,-100:100:200j]
e11=eg[0,:,:]
e22=eg[1,:,:]
e12=np.ones(np.shape(e11))

### Parameters 
e=1.3846 ## Aspect ratio of the ellipse
recip_e2=1.0/(e**2.0) #e^-2
SEAICE_strength=2.75e4 
tnsFac=0.1 # adding tensile strength, proportinnal to P*
press0=SEAICE_strength # in this case, we assume A=1 and H=1
SEAICEpressReplFac=1.0 # using the 

## Morh-Coulomb parameters
SEAICEmcMu=0.7
SEAICEmcC=SEAICEmcMu*tnsFac*SEAICE_strength

### Using smooth min and max to improve convergence
## For Delta
SEAICE_DELTA_SMOOTHREG=False
deltaMinSq=1e-20
## For Zeta
SEAICE_ZETA_SMOOTHREG=False
SEAICE_zetaMaxFac=2.5e8
ZMAX=SEAICE_zetaMaxFac*SEAICE_strength
ZMIN=1e-20
## For the MC truncating
SEAICE_MC_SMOOTHMIN=False


### Computing simpler variables
ep=e11+e22
em=e11-e22
e12Csq=e12**2.0

### Computing Delta
deltaCsq=ep*ep+recip_e2*(em*em+4.0*e12Csq)
deltaC=np.sqrt(deltaCsq)
## with a sqrt function
if SEAICE_DELTA_SMOOTHREG :
  deltaCreg = np.sqrt(deltaCsq + deltaMinSq)
## with a sharp min function
else:
  deltaCreg = np.sqrt(np.maximum(deltaCsq,deltaMinSq))


### Computing Zeta
## with a smooth tanh function
if SEAICE_ZETA_SMOOTHREG:
  argTmp = np.exp(-1./(deltaCreg*SEAICE_zetaMaxFac))
  zeta=ZMAX*(1.-argTmp)/(1.+argTmp)*(1.+tnsFac)
## with a sharp min function
else:
  zeta=press0*(1.+tnsFac)/(2.*deltaCreg)
  zeta=np.minimum(zeta,ZMAX)
zeta=np.maximum(zeta,ZMIN)

### Computing eta
eta=zeta/e**2

### Computing pressure pressure 
press=(press0*(1.-SEAICEpressReplFac)+2.*zeta*deltaC*SEAICEpressReplFac/(1.+tnsFac))*(1.-tnsFac)

### Troncated Ellipse Method 
## Computing the maximum eta to troncate ellipse
etaDen=em*em+4.*e12Csq
etaDen=np.sqrt(deltaMinSq + etaDen)
# Mu is the slope of the MC part, C is the cohesion
etaMax=(SEAICEmcMu*(0.5*press-zeta*ep)+SEAICEmcC)/etaDen 

## Troncating the ellipse 
# with smooth minimum function 
if SEAICE_MC_SMOOTHMIN :
  eta_MC=smoothmin_good(eta,etaMax,1.0)
# with a sharp min
else:
  eta_MC=np.minimum(eta,etaMax)


### Computing invariant stresses
## Ellipse
sigI=(zeta*ep-press/2.)/press0
sigII=(eta*np.sqrt(em**2.+4.*e12**2.))/press0
## Mohr-Coulomb
sigII_MC=(eta_MC*np.sqrt(em**2+4*e12**2))/press0

### Computing principal stresses
## Ellipse
sig1=sigI-sigII
sig2=sigI+sigII
## Mohr-Coulomb
sig1_MC=sigI-sigII_MC
sig2_MC=sigI+sigII_MC


### Plotting the yield curve
## Ellipse in red
## MC in blue
## theoretical ellipse in black with verticals and axis

t=tnsFac
f=1.0

fig1=plt.figure(1)
ax1 = fig1.gca()
elli=Ellipse(xy=((-f+t)/2,0), width=(f+t)/e, height=(f+t), angle=-90,edgecolor='b', fc='None', lw=0.5)
ax1.add_patch(elli)
plt.grid()
plt.plot(sigI,sigII,'.r')
plt.plot(sigI,sigII_MC,'.b')
xlim=ax1.get_xlim()
ylim=ax1.get_ylim()
ax1.set(xlim=xlim, ylim=ylim)
cordx=[min(xlim[0],ylim[0]),min(xlim[1],ylim[1])]
cordy=[min(xlim[0],ylim[0]),min(xlim[1],ylim[1])]
plt.plot(cordx,cordy,'-k',lw=1)
coradx=[-min(-xlim[0],ylim[1]),min(xlim[1],-ylim[0])]
corady=[min(-xlim[0],ylim[1]),-min(xlim[1],-ylim[0])]
plt.plot(coradx,corady,'-k',lw=1)
plt.axvline(linewidth=1, color='k')
plt.axhline(linewidth=1, color='k')
plt.axes().set_aspect('equal', 'datalim')


fig2=plt.figure(2)
ax2 = fig2.gca()
elli=Ellipse(xy=((-f+t)/2,(-f+t)/2), width=(f+t)/e*np.sqrt(2), height=(f+t)*np.sqrt(2), angle=-45,edgecolor='b', fc='None', lw=0.5)
ax2.add_patch(elli)
plt.grid()
plt.plot(sig1,sig2,'.r')
plt.plot(sig1,sig2_MC,'.b')
xlim=ax2.get_xlim()
ylim=ax2.get_ylim()
ax2.set(xlim=xlim, ylim=ylim)
cordx=[min(xlim[0],ylim[0]),min(xlim[1],ylim[1])]
cordy=[min(xlim[0],ylim[0]),min(xlim[1],ylim[1])]
plt.plot(cordx,cordy,'-k',lw=1)
coradx=[-min(-xlim[0],ylim[1]),min(xlim[1],-ylim[0])]
corady=[min(-xlim[0],ylim[1]),-min(xlim[1],-ylim[0])]
plt.plot(coradx,corady,'-k',lw=1)
plt.axvline(linewidth=1, color='k')
plt.axhline(linewidth=1, color='k')
plt.axes().set_aspect('equal', 'datalim')

plt.show()