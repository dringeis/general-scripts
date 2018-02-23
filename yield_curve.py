#!/usr/bin python
import numpy as np
from smoothabs import smoothmin_good
from matplotlib.patches import Ellipse
import pylab as plt

### Creating fake data
eg=np.mgrid[-1:1:200j,-1:1:200j]
e11=eg[0,:,:]
e22=eg[1,:,:]
e12=np.zeros(np.shape(e11))

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
SEAICE_DELTA_SMOOTHREG=True
deltaMinSq=1e-20
## For Zeta
SEAICE_ZETA_SMOOTHREG=True
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

### Computing stresses
## Computing stresses as in MITgcm
# elliptical yield curve
sig11=zeta*ep+eta*em-press/2.0
sig22=zeta*ep-eta*em-press/2.0
sig12=2*e12*eta

sigp=sig11+sig22
sigm=sig11-sig22
sigTmp=np.sqrt(sigm**2+4*sig12)

sig1=0.5*(sigp+sigTmp)*1/press0
sig2=0.5*(sigp-sigTmp)*1/press0

sigI=0.5*(sig1+sig2)
sigII=0.5*(sig1-sig2)

# MC yield curve
sig11_MC=zeta*ep+eta_MC*em-press/2.0
sig22_MC=zeta*ep-eta_MC*em-press/2.0
sig12_MC=2*e12*eta_MC

sigp_MC=sig11_MC+sig22_MC
sigm_MC=sig11_MC-sig22_MC
sigTmp_MC=np.sqrt(sigm_MC**2+4*sig12_MC)

sig1_MC=0.5*(sigp_MC+sigTmp_MC)*1/press0
sig2_MC=0.5*(sigp_MC-sigTmp_MC)*1/press0

sigI_MC=0.5*(sig1_MC+sig2_MC)
sigII_MC=0.5*(sig1_MC-sig2_MC)


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
plt.plot(sigI,sigII,'.g')
plt.plot(sigI_MC,sigII_MC,'.c')
plt.plot(sigI,-sigII,'.g')
plt.plot(sigI_MC,-sigII_MC,'.c')
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
plt.plot(sig1,sig2,'.g')
plt.plot(sig1_MC,sig2_MC,'.c')
plt.plot(sig2,sig1,'.g')
plt.plot(sig2_MC,sig1_MC,'.c')
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