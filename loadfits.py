import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from scipy import optimize

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
                print(year, band, " matched " filename)
                hdul = fits.open(root + filename)
                images.append(hdul[0].data)
                seeing_pixels.append(hdul[0].header["L1SEEING"])
                seeing_arcsec.append(hdul[0].header["L1SEESEC"])
    return(images)

def gaussian(height, center_x, center_y, width_x, width_y):
    """
    Returns a gaussian function with the given parameters.
    """
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*np.exp(-(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

def moments(data):
    """
    Returns (height, x, y, width_x, width_y) the gaussian parameters of a 2D
    distribution by calculating its moments.
    """
    total = data.sum()
    X, Y = np.indices(data.shape)
    x = (X * data).sum() / total
    y = (Y * data).sum() / total
    col = data[:,int(y)]
    width_x = np.sqrt(np.abs((np.arange(col.size) - y)**2 * col).sum() / col.sum())
    row = data[int(x), :]
    width_y = np.sqrt(np.abs((np.arange(row.size) - x)**2 * row).sum() / row.sum())
    height = data.max()
    return height, x, y, width_x, width_y

def fitgaussian(data):
    """
    Returns (height, x, y, width_x, width_y) the gaussian parameters of a 2D
    distribution found by a fit.
    """
    params = moments(data)
    errorfunction = lambda p: np.ravel(gaussian(*p)(*np.indices(data.shape)) - data)
    p, success = optimize.leastsq(errorfunction, params)
    return(p)

def radians(coordinates):
    """
    Receives a tuple of strings, containing RA and DEC coordinates in the format:
    ("NN:NN:NN.NN","+NN:NN:NN.NN")
    And returns a tuple of the equivalent values in radians.
    """
    ra_str, dec_str = coordinates[0].split(":"), coordinates[1].split(":")
    ra = np.pi * (ra_str[0] / 12 + ra_str[1] / 720 + ra_str[2] / 43200)
    dec = np.pi * (dec_str[0] / 12 + dec_str[1] / 720 + dec_str[2] / 43200)
    return((ra, dec))

def degrees(coordinates):
    """
    Receives a tuple of strings, containing RA and DEC coordinates in the format:
    ("NN:NN:NN.NN","+NN:NN:NN.NN")
    And returns a tuple of the equivalent values in degrees.
    """
    ra_str, dec_str = coordinates[0].split(":"), coordinates[1].split(":")
    ra = 360 * (ra_str[0] / 24 + ra_str[1] / 1440 + ra_str[3] / 86400)
    dec = 360 * (dec_str[0] / 24 + dec_str[1] / 1440 + dec_str[3] / 86400)
    return((ra, dec))

def bad_centroid(cutout):
    """
    Receives a cutout of an image, and returns a tuple of the coordinates of the
    brightest pixel in that image, in pixels.
    """
    x_sum = np.sum(cutout, axis=0)
    y_sum = np.sum(cutout, axis=1)
    x_max = np.where(x_sum == max(x_sum))[0][0]
    y_max = np.where(y_sum == max(y_sum))[0][0]
    return((x_max, y_max))

def bad_centroid_2(cutout, **kwargs):
    """
    Recieves an argument of type ndarray and returns a tuple of the weighted
    mean centroid of the object contained in the cutout, omitting the values
    which fall beneath half the maximum of the image.
    """
    x_sum = np.sum(cutout, axis=0)
    y_sum = np.sum(cutout, axis=1)
    x_fwh = np.where(x_sum >= 0.5 * max(x_sum))
    y_fwh = np.where(y_sum >= 0.5 * max(y_sum))
    x_avg = np.average(x_fwh[0], weights=x_sum[x_fwh])
    y_avg = np.average(y_fwh[0], weights=y_sum[y_fwh])
    plt.imshow(cutout, cmap='viridis', origin='lower')
    plt.scatter(x_avg, y_avg, s=2, c='red', marker='o')
    plt.show()
    return((int(np.floor(x_avg)), int(np.floor(y_avg))))

def annuli(array, center, radii):
    """
    Receives and array and returns a tuple of three masked annuli from the
    input array about a given center.
    """
    return inner_annulus, middle_annulus, outer_annulus

def hybrid_centroid(image, object_coords, object_size, **kwargs):
    """
    Receives a FITS image, proper coordinates of an object in the image, the size
    """

    obj_coords = degrees(object_coords)
    ref_coords = degrees((image[0].header["RA"], image[0].header["DEC"]))

    deg_offset = obj_coords - ref_coords

    ra_pix = image[0].header["CDELT1"]
    dec_pix = image[0].header["CDELT2"]
    pix_offset = (ra_pix * deg_offset[0], dec_pix * deg_offset[1])

    obj_guess = (pix_offset) + (image[0].header["CRPIX1"], image[0].header["CRPIX2"])

    image_data = image[0].data

    annuli = annulii(image_data, obj_guess, obj_size)

    x_sum = np.sum(annulii[0], axis=0)
    y_sum = np.sum(annulii[0], axis=1)
    # plt.plot(range(len(x_sum)), x_sum)
    # plt.show()
    # plt.plot(range(len(y_sum)), y_sum)
    # plt.show()
    x_avg = np.average(range(x_sum.size), weights=x_sum)
    y_avg = np.average(range(y_sum.size), weights=y_sum)

    return((x_avg, y_avg))

def weighted_mean_2D(cutout,**kwargs):
    """
        Recieves an argument of type ndarray and returns a tuple of the weighted
        mean centroid of the object contained in the cutout.
    """
    x_sum = np.sum(cutout, axis=0)
    y_sum = np.sum(cutout, axis=1)
    # plt.plot(range(len(x_sum)), x_sum)
    # plt.show()
    # plt.plot(range(len(y_sum)), y_sum)
    # plt.show()
    x_avg = np.average(range(x_sum.size), weights=x_sum)
    y_avg = np.average(range(y_sum.size), weights=y_sum)
    plt.imshow(cutout, cmap='viridis', origin='lower')
    plt.scatter(x_avg, y_avg, s=2, c='red', marker='o')
    plt.show()
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
    centroid = kwargs.get("centroid")
    # Get lists of all the x and y centroids.
    x_centroids, y_centroids = [], []
    for image in image_stack:
        x_centroids.append(centroid(image[x:x+dx,y:y+dy],floor=True)[0])
        y_centroids.append(centroid(image[x:x+dx,y:y+dy],floor=True)[1])
    x_ref, y_ref = min(x_centroids), min(y_centroids)
    x_max_offset = max(x_centroids) - min(x_centroids)
    y_max_offset = max(y_centroids) - min(y_centroids)
    # Create new list of image arrays with offset.
    aligned_image_size = ()
    aligned_image_stack = []
    for image in image_stack:
        aligned_image = np.zeros((image.shape[0]+x_max_offset, image.shape[1]+y_max_offset))
        x_image_offset = centroid(image[x:x+dx,y:y+dy],floor=True)[0] - x_ref
        y_image_offset = centroid(image[x:x+dx,y:y+dy],floor=True)[1] - y_ref
        print(x_image_offset, y_image_offset)
        aligned_image[x_image_offset:x_image_offset+image.shape[0],y_image_offset:y_image_offset+image.shape[1]] = image
        aligned_image_stack.append(aligned_image)
    return(aligned_image_stack)

def stack(aligned_image_stack):
    # plt.plot(range(len(x_sum)), x_sum)
    # plt.show()
    # plt.plot(range(len(y_sum)), y_sum)
    # plt.show()
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
    ax_array[0].imshow(image_r[400:600,400:600], cmap=cmap, origin='lower', norm=LogNorm())
    ax_array[1].imshow(image_g[400:600,400:600], cmap=cmap, origin='lower', norm=LogNorm())
    ax_array[1].scatter(524-400, 503-400, s=2, c='red', marker='o')
    ax_array[1].scatter(472-400, 496-400, s=2, c='red', marker='o')
    ax_array[1].scatter(487-400, 473-400, s=2, c='red', marker='o')
    ax_array[2].imshow(image_u[400:600,400:600], cmap=cmap, origin='lower', norm=LogNorm())
    plt.show()

def rgb(image_r, image_g, image_b):
    """
        Recieves three arrays of equal size. Maps these values to RGB values
        using the Lupton algorithm and displays the resulting image.
        # TODO: Retrieve the source for this algorithm.
    """
    from astropy.visualization import make_lupton_rgb
    rgb_image = make_lupton_rgb(image_r, image_g, image_b, Q=10, stretch=1000.)
    plt.imshow(rgb_image)
    plt.show()

def hist(list1, list2):
    """
    Recieves two lists of seeing values, creates a bin histogram of each and
    displays tthe resulting plots.
    """
    figure, ax_array = plt.subplots(1,2)
    ax_array[0].hist(list1, range=(3,10))
    ax_array[1].hist(list2, range=(2,2.6))
    plt.show()

seeing_pixels, seeing_arcsec = [] , []

def main():
    # Empty list for collecting stacked images in three bands.
    rgu_images = []
    # Define the cutout region containing the reference object.
    x, y, dx, dy = 450, 450, 75, 60
    # Define the proper location of the object, for alignment.
    # proper_coords = radians(("09:45:11.08","17:45:44.80"))

    for  year in range(2012,2018):
        for band in ["R", "G", "U"]:
            unaligned_images = load_fits(path="data/SDSSJ094511-P1-images/", year=str(year), band=band)
            aligned_images = align(unaligned_images, cutout=(x,y,dx,dy), centroid=bad_centroid_2)
            stacked_image = stack(aligned_images)
            rgu_images.append(stacked_image)

        plot(rgu_images[0], rgu_images[1], rgu_images[2], 'viridis')
        rgb(rgu_images[0][:1024,:1024], rgu_images[1][:1024,:1024], rgu_images[2][:1024,:1024])
        hist(seeing_pixels, seeing_arcsec)
main()
