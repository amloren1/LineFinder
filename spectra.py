import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import scipy
from scipy import interp
from scipy import signal
import math
import pdb
import sys




readpath_lb1 = r'data/lb1v15010_x1dsum.fits' #object 1
#readpath_lcd = 'data/lcd905010_x1dsum.fits'  #objext 2



obj2= 'lcd905010_x1dsum.fits'
obj1='lb1v15010_x1dsum.fits'

a = fits.open(readpath_lb1)

a.info()

hdu = a[1]

#with open('test.txt', 'w') as f: f.write(repr(hdu.header))





w = a[1].data['WAVELENGTH']
f = a[1].data['FLUX']



print w.shape
print f.shape







def moving_average(a, n=8) :
    #n is the number of points to average
    nn = len(a)/n
    print 'from {} to {}'.format(len(a),nn) 
    csum = np.cumsum(a, dtype=float)
    print csum
    ret = np.zeros(nn)
    lo =0
    for i in range(nn):
        hi = i*n+(n-1)
        #print lo
        #print hi
        ret[i] = (csum[hi]-csum[lo])/n
        lo =hi
        #print ret 
    return ret

test = np.arange(0,12)


ww = moving_average(w[0])
ff = moving_average(f[0])


def gaussfold(lam, flux, fwhm):
    #more reading info http://www.igsfddsaforexchange.com/node/2985
    lammin=min(lam)
    lammax=max(lam)
    dlambda=fwhm/17.
    interlam=np.arange(lammin,lammax,dlambda)
    interflux=interp(interlam,lam,flux)
    fwhm_pix=fwhm/dlambda
    window=math.floor(17*fwhm_pix)

    #constructing a normalized gaussian who's width is FWHM
    std=0.5*fwhm_pix/math.sqrt(2*math.log(2))
    gauss1=scipy.signal.gaussian(window,std=std)
    gauss=gauss1/np.trapz(gauss1)

    # pdb.set_trace()
    #convovle flux array with gaussian
    fold=np.convolve(interflux,gauss,mode='same')
    #interpolate back to original grid
    fluxfold=interp(lam,interlam,fold)
    #pdb.set_trace()
    return fluxfold

#this running/boxvar mean thing doesn't work so well here
# wbin = moving_average(w[0],100)
# fbin = []
# fbin.append(np.average(f[0][np.where((w[0] >= w[0][0]) & (w[0] < wbin[0]))]))
# for i in range(len(wbin[1:])):
#     loc1 = np.where((w[0] >= wbin[0+i]) & (w[0] < wbin[1+i]))
#     fbin.append(np.average(f[0][loc1]))

#R is the new spectral resolution, R = lambda / delta lamda
R = 20000
fwhm = 1./R * np.mean(w[0])
fbin = gaussfold(w[0],f[0],fwhm)

#


data = np.vstack((ww,ff))
data = np.transpose(data)
np.savetxt('lb1v15010_test.dat', data, delimiter=' ', comments='#', header = 'WAVELENGTH FLUX')


#sys.exit()

##################
#plots
##################

fig, (ax1,ax2) = plt.subplots(2,1)

ax1.plot(ww,ff,color='b', label = 'data')
ax2.plot(w[0],fbin,color='k', label='reduction')
ax2.set_xlabel('Wavelength')
ax1.set_ylabel('Flux')
ax2.set_ylabel('Flux')
# plot(w[1],f[1],color='b')
ax1.set_xlim([1405,1445])
ax2.set_xlim([1405,1445])
ax1.legend()
ax2.legend()
plt.show()
