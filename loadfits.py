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

def plot():

    from matplotlib.colors import LogNorm

    figure, ax_array = plt.subplots(2,3)
    ax_array[0,0].imshow(np.rot90(image_2012_r, k=1)[475:525,475:525], cmap='viridis', origin='lower', norm=LogNorm())
    ax_array[0,1].imshow(np.rot90(image_2012_r, k=1)[475:525,475:525], cmap='viridis', origin='lower', norm=LogNorm())
    ax_array[0,2].imshow(np.rot90(image_2012_r, k=1)[475:525,475:525], cmap='viridis', origin='lower', norm=LogNorm())
    ax_array[1,0].imshow(np.rot90(image_2012_r, k=1)[475:525,475:525], cmap='viridis', origin='lower', norm=LogNorm())
    ax_array[1,1].imshow(np.rot90(image_2012_r, k=1)[475:525,475:525], cmap='viridis', origin='lower', norm=LogNorm())
    ax_array[1,2].imshow(np.rot90(image_2012_r, k=1)[475:525,475:525], cmap='viridis', origin='lower', norm=LogNorm())
    plt.show()

def rgb(image_r, image_g, image_b):
    """
        Recieves three arrays of equal size. Maps these values to RGB values
        using the Lupton algorithm and displays the resulting image.
        # TODO: Retrieve the source for this algorithm.
    """
    from astropy.visualization import make_lupton_rgb
    rgb_image = make_lupton_rgb(image_r, image_g, image_b, Q=10, stretch=1000.)
    plt.imshow(image)
    plt.show()

def weighted_mean_2D(cutout,**kwargs):
    """
        Recieves an argument of type ndarray and returns a tuple of the weighted
        mean centroid of the object contained in the cutout.
    """
    x_sum = np.sum(cutout, axis=0)
    y_sum = np.sum(cutout, axis=1)
    x_avg = np.average(range(x_sum.size), weights=x_sum)
    y_avg = np.average(range(y_sum.size), weights=y_sum)
    if kwargs.get("floor") == True:
        return((np.floor(x_avg), np.floor(y_avg)))
    else:
        return((x_avg, y_avg))

def align(image_stack, **kwargs):
    """
    Recieves a list of image arrays and some "cutout" range containing a common
    objecct to use for alignment of the image stack.

    Returns a list of image arrays of different size, aligned, and with zero
    borders where the image has been shifted.
    """

    cutout_range = kwargs.get("cutout")

    # Get lists of all the x and y centroids.
    for image in image_stack:
        x_centroids.append(weighted_mean_2D(image[cutout_range])[0]):
        y_centroids.append(weighted_mean_2D(image[cutout_range])[1]):
    x_ref, y_ref = x_centroids.min(), y_centroids.min()
    x_max, y_max = x_centroids.max(), y_centroids.max()

    # Create new list of image arrays with offset.
    aligned_image_stack = []
    for image in image_stack:
        aligned_image = np.zeros(image.size[0]+x_max, image.size[1]+y_max)
    return(aligned_image_stack)

def main():

    image_2012_r = fits.open("data/fits/20120312_39_R100.fits")[0].data[475:525,475:525]
    image_2012_g = fits.open("data/fits/20120312_38_G100.fits")[0].data[475:525,475:525]
    image_2012_b = fits.open("data/fits/20120312_40_U300.fits")[0].data[475:525,475:525]
    image_2017_r = fits.open("data/fits/20170314_20_R.fits")[0].data[475:525,475:525]
    image_2017_g = fits.open("data/fits/20170314_19_G.fits")[0].data[475:525,475:525]
    image_2017_b = fits.open("data/fits/20170314_21_U.fits")[0].data[475:525,475:525]

    images_2012 = [image_2012_r,image_2012_g,image_2012_b]
    images_2017 = [image_2017_r,image_2017_g,image_2017_b]

    # plt.imshow(np.rot90(image_stack, k=1)[495:505,495:505], cmap='viridis', origin='lower', norm=LogNorm())

main()
