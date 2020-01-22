import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np

"""
CCDXIMSI=                 1024 / [pixels] Imaging pixels
CCDYIMSI=                 1024 / [pixels] Imaging pixels
CCDXBIN =                    2 / [pixels] X binning factor
CCDYBIN =                    2 / [pixels] Y binning factor
CCDXPIXE=            0.0000135 / [m] Size of pixels, in X:13.5um
CCDYPIXE=            0.0000135 / [m] Size of pixels, in Y:13.5um
CCDSCALE=                0.279 / [arcsec/binned pixel] Scale of binned image on
RA      = '09:45:10.866'       / World coordinate at the reference pixel
DEC     = '+17:45:48.17'       / World coordinate at the reference pixel
RADECSYS= 'FK5     '           / [FK4, FK5] Fundamental coord system
L1SEEING=             4.350539 / [pixels] frame seeing in pixels
L1SEESEC=             1.213801 / [arcsec] frame seeing in arcsec
SEEING  =             4.350539 / [pixels] frame seeing in pixels
CTYPE1  = 'RA---TAN'           / WCS projection
CTYPE2  = 'DEC--TAN'           / WCS projection
CRPIX1  =                 512. / WCS reference pixel coordinate
CRPIX2  =                 512. / WCS reference pixel coordinate
CRVAL1  =        146.295274426 / [degrees] World coordinate at the ref pix
CRVAL2  =          17.76338016 / [degrees] World coordinate at the ref pix
WRA     = '9:45:11.075'        / Original RA value from TCS before WCS fit
WDEC    = '+17:45:44.89'       / Original DEC value from TCS before WCS fit
WRADECSY= 'FK5     '           / [FK4,FK5] Original Fundamental system before WC
WROTSKY =             270.0001 / [deg] Original sky PA from TCS before WCS fit
EPOCH   =                2000. / Epoch of coordinate
CDELT1  =       -7.7508973E-05 / [degrees/pixel]
CDELT2  =        7.7508973E-05 / [degrees/pixel]
CROTA1  =            90.361783 / [degrees]
CROTA2  =            90.361783 / [degrees]
"""

def weighted_mean_2D(array):
    """
        Recieves an argument of type ndarray and returns a tuple of the weighted
        mean centroid of the object.
    """
    x_sum = np.sum(array, axis=0)
    y_sum = np.sum(array, axis=1)
    print(x_sum)
    print(y_sum)
    x_avg = np.average(range(x_sum.size), weights=x_sum)
    y_avg = np.average(range(y_sum.size), weights=y_sum)
    # x_avg = 0.0
    # for i in range(len(x_sum)):
    #     x_avg += float(i) * x_sum[i] / np.sum(x_sum)
    # y_avg = 0.0
    # for j in range(len(y_sum)):
    #     y_avg += float(j) * x_sum[j] / np.sum(y_sum)

    return((x_avg, y_avg))

image_2012_r = fits.open("data/fits/20120312_39_R100.fits")[0].data
image_2012_g = fits.open("data/fits/20120312_38_G100.fits")[0].data
image_2012_b = fits.open("data/fits/20120312_40_U300.fits")[0].data
image_2017_r = fits.open("data/fits/20170314_20_R.fits")[0].data
image_2017_g = fits.open("data/fits/20170314_19_G.fits")[0].data
image_2017_b = fits.open("data/fits/20170314_21_U.fits")[0].data

# hdul = fits.open("data/fits/20120312_38_G100.fits")
# hdul.info()

#subplot(r,c) provide the no. of rows and columns

from matplotlib.colors import LogNorm

figure, ax_array = plt.subplots(2,3)
ax_array[0,0].imshow(np.rot90(image_2012_r, k=1), cmap='viridis', origin='lower', norm=LogNorm())
ax_array[0,1].imshow(np.rot90(image_2012_r, k=1), cmap='viridis', origin='lower', norm=LogNorm())
ax_array[0,2].imshow(np.rot90(image_2012_r, k=1), cmap='viridis', origin='lower', norm=LogNorm())
ax_array[1,0].imshow(np.rot90(image_2012_r, k=1), cmap='viridis', origin='lower', norm=LogNorm())
ax_array[1,1].imshow(np.rot90(image_2012_r, k=1), cmap='viridis', origin='lower', norm=LogNorm())
ax_array[1,2].imshow(np.rot90(image_2012_r, k=1), cmap='viridis', origin='lower', norm=LogNorm())
plt.show()

figure, ax_array = plt.subplots(2,3)
ax_array[0,0].imshow(np.rot90(image_2012_r, k=1)[475:525,475:525], cmap='viridis', origin='lower', norm=LogNorm())
ax_array[0,1].imshow(np.rot90(image_2012_r, k=1)[475:525,475:525], cmap='viridis', origin='lower', norm=LogNorm())
ax_array[0,2].imshow(np.rot90(image_2012_r, k=1)[475:525,475:525], cmap='viridis', origin='lower', norm=LogNorm())
ax_array[1,0].imshow(np.rot90(image_2012_r, k=1)[475:525,475:525], cmap='viridis', origin='lower', norm=LogNorm())
ax_array[1,1].imshow(np.rot90(image_2012_r, k=1)[475:525,475:525], cmap='viridis', origin='lower', norm=LogNorm())
ax_array[1,2].imshow(np.rot90(image_2012_r, k=1)[475:525,475:525], cmap='viridis', origin='lower', norm=LogNorm())
plt.show()

image_stack = image_2012_r + image_2012_g + image_2012_b

print(weighted_mean_2D(image_stack))

plt.imshow(np.rot90(image_stack, k=1)[495:505,495:505], cmap='viridis', origin='lower', norm=LogNorm())
plt.colorbar()
plt.show()
