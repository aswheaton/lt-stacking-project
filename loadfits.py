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

def load_fits(**kwargs):
    """
    Receives a directory path and .fits filename parameters. Parses the
    directory for files matching the naming parameters and loads matched files
    into a list of astropy.fits objects. Returns the list.
    """
    path = kwargs.get("path")
    year = kwargs.get("year")
    band = kwargs.get("band")

    import os

    images = []
    for root, dirs, files in os.walk(path):
        for filename in files:
            if year in filename and band in filename:
                images.append(fits.open(root + filename)[0].data)
    return(images)

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
        return((int(np.floor(x_avg)), int(np.floor(y_avg))))
    else:
        return((x_avg, y_avg))

def align(image_stack, **kwargs):
    """
    Recieves a list of image arrays and some "cutout" range containing a common
    object to use for alignment of the image stack.

    Returns a list of image arrays of different size, aligned, and with zero
    borders where the image has been shifted.
    """
    x, y, dx, dy = kwargs.get("cutout")
    # Get lists of all the x and y centroids.
    x_centroids, y_centroids = [], []
    for image in image_stack:
        x_centroids.append(weighted_mean_2D(image[x:x+dx,y:y+dy],floor=True)[0])
        y_centroids.append(weighted_mean_2D(image[x:x+dx,y:y+dy],floor=True)[1])
    x_ref, y_ref = min(x_centroids), min(y_centroids)
    x_max_offset = max(x_centroids) - min(x_centroids)
    y_max_offset = max(y_centroids) - min(y_centroids)
    # Create new list of image arrays with offset.
    aligned_image_size = ()
    aligned_image_stack = []
    for image in image_stack:
        aligned_image = np.zeros((image.shape[0]+x_max_offset, image.shape[1]+y_max_offset))
        x_image_offset = weighted_mean_2D(image[x:x+dx,y:y+dy],floor=True)[0] - x_ref
        y_image_offset = weighted_mean_2D(image[x:x+dx,y:y+dy],floor=True)[1] - y_ref
        aligned_image[x_image_offset:x_image_offset+image.shape[0],y_image_offset:y_image_offset+image.shape[1]] = image
        aligned_image_stack.append(aligned_image)
    return(aligned_image_stack)

def stack(aligned_image_stack):
    """
        Receives a list of aligned images and returns their sum.
    """
    # Check that the aligned images to be stacked have matching dimensions.
    for image in aligned_image_stack:
        if image.shape != aligned_image_stack[0].shape:
            print("Aligned image dimensions do not match!")
            break
    # If all dimensions match, initialise an empty array with those dimensions
    # into which aligned images are stacked.
    stacked_image = np.zeros(aligned_image_stack[0].shape)
    for image in aligned_image_stack:
        stacked_image += image
    return(stacked_image)

def plot(image_r, image_g, image_u, cmap):
    """
        Recieves image arrays in three bands and plots them according to a given
        colour map, on a log scale.
    """
    from matplotlib.colors import LogNorm
    figure, ax_array = plt.subplots(1,3)
    ax_array[0,0].imshow(image_r, cmap=cmap, origin='lower', norm=LogNorm())
    ax_array[0,1].imshow(image_g, cmap=cmap, origin='lower', norm=LogNorm())
    ax_array[0,2].imshow(image_u, cmap=cmap, origin='lower', norm=LogNorm())
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

def main():

    unaligned_images = loadfits(path="data/SDSSJ094511-P1-images/", year="2012", band="R")
    aligned_images = align(unaligned_images, cutout=[475:525,475:525])
    stacked_image = stack(aligned_images)


    # plt.imshow(np.rot90(image_stack, k=1)[495:505,495:505], cmap='viridis', origin='lower', norm=LogNorm())

main()
