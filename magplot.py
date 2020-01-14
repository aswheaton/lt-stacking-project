# J094511 light curve

#%% set up
import numpy as np
from matplotlib import pyplot as plt
from math import *
from astropy.io import ascii 

#%% info parameters
name='J094511'
z=0.758

## plot params
tmin=52000  # MJD
tmax=58500
magmin=23.0  # g or r mag as appropriate
magmax=19.0

Voff=0.3  # offset for CRTS to g-band

#%% read in data from txt files
LT=ascii.read("lc4-LT.txt")  
tLT=LT['MJD']
gLT=LT['g']
rLT=LT['r']
gLTerr=LT['gerr']
rLTerr=LT['rerr']

dvog=ascii.read('lc4-dvo-g.txt')  
tdvog=dvog['MJD']
gdvo=dvog['mag']
gdvoerr=dvog['magerr']

dvor=ascii.read('lc4-dvo-r.txt')  
tdvor=dvor['MJD']
rdvo=dvor['mag']
rdvoerr=dvor['magerr']

crts=ascii.read('lc4-crts.txt')  
tcrts=crts['MJD']
Vcrts=crts['V']
gcrts=Vcrts+Voff
gcrtserr=crts['Verr'] # note use Verr as gerr

sdss=ascii.read('lc4-sdr9.txt')  
tsdss=sdss['MJD']
gsdss=sdss['g']
gsdsserr=sdss['gerr']
rsdss=sdss['r']
rsdsserr=sdss['rerr']

spec=ascii.read('lc4-spec.txt')
tspec=spec['MJD']

#%% make plot
fig1=plt.figure()
ax1=fig1.gca()
ax1.axis([tmin,tmax,magmin,magmax])
ax1.set_title(name+' z='+str(z), size=16)
ax1.set_xlabel('MJD', size=16)
ax1.set_ylabel('gmag', size=16)
#ax1.set_ylabel('rmag', size=16)
ax1.tick_params(labelsize=16)

## g-band version
ax1.errorbar(tLT,gLT,fmt='o',yerr=gLTerr)
ax1.errorbar(tdvog,gdvo,fmt='o',yerr=gdvoerr)
ax1.errorbar(tsdss,gsdss,fmt='o',yerr=gsdsserr)

## r-rand version
#ax1.errorbar(tLT,rLT,fmt='o',yerr=rLTerr)
#ax1.errorbar(tdvor,rdvo,fmt='o',yerr=rdvoerr)
#ax1.errorbar(tsdss,rsdss,fmt='o',yerr=rsdsserr)

## CRTS, labelled g but can shift to R by changing Voff
ax1.errorbar(tcrts,gcrts,fmt='o',yerr=gcrtserr)

# label spectral epochs
for i in range(0,len(tspec)):
    ax1.plot([tspec[i],tspec[i]],[magmin,magmax],'k--')


plt.show()






