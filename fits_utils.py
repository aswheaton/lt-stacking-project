#! /usr/bin/env python

"""
fits-utils module

This module contains various functions for the reduction, alignment and stacking
of raw astronomical frames using the FITS file formatting standard. It is an
extension of the astropy module, and utilises many of its basic functions to
perform more complicated operations.

The raw files are stored in the "dat/" folder, which is read-only; science
frames are written to and read from in the "sci/" folder; temporary files
(including master dark and flat fields) are written to and read from the "tmp/"
directory; and stacked images are written to and read from the "sta/" directory.

This module utilises the configparser module in order to generate a "config.ini"
file. This file can be populated by the user with a target ID, standard star ID,
and observing bands. The module can then use these parameters to parse the
"raw/" directory and sort FITS files into lists with specific integration times
and bands. The images can then be reduced and written out to FITS files.

By convention, imported:

import fits-utils as fu
"""

import numpy as np
import matplotlib.pyplot as plt
import configparser

from matplotlib.colors import LogNorm
from scipy.stats import mode
from scipy.ndimage import gaussian_filter
from pathlib import Path
from astropy.io import fits
from os import walk

# Force fonts to be nice.
# from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern Sans serif']})
# rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
# rc('text', usetex=True)

def load_fits(**kwargs):
    """
    Receives a directory path and .fits filename parameters. Parses the
    directory for files matching the naming parameters and loads matched files
    into a list of astropy.fits objects. Returns the list.
    """
    path = kwargs.get("path")
    year = kwargs.get("year")
    band = kwargs.get("band")
    images = []
    for root, dir, files in walk(path):
        for filename in files:
            if year in filename and band in filename:
                print(year, band, " matched ", filename)
                with fits.open(root+filename) as hdul:
                    new_image = {}
                    new_image["filename"] = filename
                    new_image["data"] = hdul[0].data
                    new_image["year"] = year
                    new_image["band"] = band
                    new_image["rot1"] = hdul[0].header["CROTA1"]
                    new_image["rot2"] = hdul[0].header["CROTA1"]

                    images.append(new_image)
    return(images)

def write_out_fits(image, filename):
    """
    Creates a header for an ndarray of reduced data and then creates a new fits
    file of this data.

    Args:
        image (ndarray): reduced data to be written to fits file.
        filename (string): name (and location) of new fits file.
    """
    hdul = fits.HDUList([fits.PrimaryHDU(image)])
    hdul.writeto(filename, overwrite=True)

def write_out_fits_2(image, filename):
    """
    Creates a header for an ndarray of reduced data and then creates a new fits
    file of this data.

    Args:
        image (ndarray): reduced data to be written to fits file.
        filename (string): name (and location) of new fits file.
    """
    hdul = fits.HDUList([fits.PrimaryHDU(image["data"])])
    hdul.writeto(filename, overwrite=True)

def get_click_coord(array):

    def onclick(click):
        global point
        # Reversed because (row, column) = (y, x)
        point = (click.ydata,click.xdata)
        return(point)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.imshow(array[400:600,400:600], cmap='viridis', origin='lower', norm=LogNorm())
    plt.show()

    floored_point = (int(np.trunc(point)[0])+400, int(np.trunc(point)[1])+400)

    return(floored_point)

def weighted_mean_2D(cutout,**kwargs):
    """
    Recieves an argument of type ndarray and returns a tuple of the weighted
    mean centroid of the object contained in the cutout.

    Args:
        cutout (ndarray): portion of fits full ndarray.
    Returns:
        x_avg (int): weighted mean of x values.
        y_avg (int): weighted mean of y values.
    """
    x_sum = np.sum(cutout, axis=0)
    y_sum = np.sum(cutout, axis=1)
    x_avg = np.average(range(x_sum.size), weights=x_sum)
    y_avg = np.average(range(y_sum.size), weights=y_sum)
    if kwargs.get("floor") == True:
        return((int(np.floor(x_avg)), int(np.floor(y_avg))))
    else:
        return((x_avg, y_avg))

def align(images, **kwargs):
    """
    Recieves a list of image arrays containing a common object to use for
    alignment of the image stack. Returns a list of image arrays of different
    size, aligned, and with zero borders where the image has been shifted.

    Args:
        images (list of dict): Frames to be aligned.

    Returns:
        aligned_images (list of dict): new frames that have been aligned and can
            be stacked.
    """
    # Which centroiding function to use.
    centroid_func = kwargs.get("centroid")
    # Boolean, whether or not to mask images for hot pixels on the detector.
    filter = kwargs.get("filter")
    # Find the centroid of the reference star in each image.
    x_centroids, y_centroids = [], []
    print("---Beginning Alignment---")
    counter = 0
    for image in images:
        counter += 1
        print("---Finding Centre {} of {}".format(counter, len(images)), end="\r")
        x_guess, y_guess = get_click_coord(image["data"])
        cutout = image["data"][x_guess-5:x_guess+5,
                               y_guess-5:y_guess+5
                               ]
        # Get the mean weighted average of the smaller cutout.
        x_new, y_new = weighted_mean_2D(cutout, floor=True)
        # Map the coordinates back to the whole image.
        x_new = x_guess + x_new - 5
        y_new = y_guess + y_new - 5
        x_centroids.append(x_new)
        y_centroids.append(y_new)
        image["XCENT"] = x_new
        image["YCENT"] = y_new

    print()
    max_pos_x, max_pos_y = max(x_centroids), max(y_centroids)
    min_pos_x, min_pos_y = min(x_centroids), min(y_centroids)
    max_dif_x, max_dif_y = max_pos_x - min_pos_x, max_pos_y - min_pos_y
    # Create new stack of aligned images using the centroid in each frame.
    aligned_images = []
    counter = 0
    for image in images:
        counter += 1
        print("---Aligning Image {} of {}".format(counter, len(images)), end="\r")
        # Determine region in which to cast the image.
        disp_x, disp_y = max_pos_x - image["XCENT"], max_pos_y - image["YCENT"]
        print("\nShifting image {} by {}, {}".format(counter,disp_x,disp_y))
        imsize_x, imsize_y = image["data"].shape[0], image["data"].shape[1]
        # Create new array containing aligned image data.
        aligned_image_data = np.zeros((imsize_x+max_dif_x, imsize_y+max_dif_y), dtype=int)
        aligned_image_data[disp_x:disp_x+imsize_x,disp_y:disp_y+imsize_y] = image["data"]
        # Create new dictionary value of frames per pixel column (one everywhere for unstacked image).
        frame_count = np.zeros((imsize_x+max_dif_x, imsize_y+max_dif_y), dtype=int)
        frame_count[disp_x:disp_x+imsize_x,disp_y:disp_y+imsize_y] = 1
        # Create new image dictionary and copy over header data from image.
        aligned_image = {"filename"    : image["filename"],
                         "data"        : aligned_image_data,
                         "frame_count" : frame_count
                         }
        # Add the new aligned image dictionary to a list to be returned.
        aligned_images.append(aligned_image)
    print("---Alignment Complete---")
    return(aligned_images)

def stack(aligned_image_stack, **kwargs):
    """
    Receives a list of aligned images and returns their summation along the axis
    of the list.

    Args:
        aligned_image_stack (list of dict): aligned frames ready to be stacked.
    Returns:
        stacked_image (dict): new combined single frame.
    """

    print("---Stacking Images---")

    # Check that the aligned images to be stacked have matching dimensions.
    for image in aligned_image_stack:
        if image["data"].shape != aligned_image_stack[0]["data"].shape:
            print("Aligned image dimensions do not match!")
            break

    # Initialise an empty array into which aligned images are stacked.
    stacked_image_data = np.zeros(aligned_image_stack[0]["data"].shape)
    counter = 0

    for image in aligned_image_stack:
        counter += 1
        print("---Stacking Image {} of {}".format(counter, len(aligned_image_stack)), end="\r")
        stacked_image_data += image["data"]
    stacked_image = {"data" : stacked_image_data}

    return(stacked_image)

def manual_centroid(image_data, **kwargs):
    plt.imshow(image_data, cmap='viridis', origin='lower', norm=LogNorm())
    plt.show()
    rotations = int(input("\nRotate cc: "))
    image_data = np.rot90(image_data, rotations)
    x_guess, y_guess = get_click_coord(image_data)
    cutout = np.array(image_data[x_guess-5:x_guess+5,
                                 y_guess-5:y_guess+5
                                 ])
    # Get the mean weighted average of the smaller cutout.
    x_new, y_new = weighted_mean_2D(cutout, floor=True)
    return((x_new, y_new))

from scipy import optimize

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
    ra = 360 * (float(ra_str[0]) / 24 + float(ra_str[1]) / 1440 + float(ra_str[2]) / 86400)
    dec = 360 * (float(dec_str[0]) / 24 + float(dec_str[1]) / 1440 + float(dec_str[2]) / 86400)
    return((ra, dec))

def wcs_offset(proper_coords, image):
    central_image_coords = (image[0].header["CRPIX1"], image[0].header["CRPIX2"])
    ra_offset = degrees(proper_coords)[0] - degrees(central_image_coords)[0]
    dec_offset = degrees(proper_coords)[1] - degrees(central_image_coords)[1]
    ra_pix = image[0].header["CDELT1"]
    dec_pix = image[0].header["CDELT2"]
    pix_offset = (ra_pix * ra_offset, dec_pix * dec_offset)
    obj_guess = (pix_offset[0] + central_image_coords[0], pix_offset[1] + central_image_coords[1])
    return(obj_guess)

def annuli_mask(array, center, radii):
    """
    Receives and array and returns a tuple of three masked annuli from the
    input array about a given center.
    """

    inner_annulus = np.zeros(6*radii, 6*radii)
    middle_annulus = np.zeros(6*radii, 6*radii)
    outer_annulus = np.zeros(6*radii, 6*radii)

    for x in range(6*radii):
        for y in range(6*radii):
            sumsqu = (x-(3*radii))**2 + (y-(3*radii))**2
            if sumsqu < radii:
                inner_annulus[x,y] = 0
                middle_annulus[x,y] = 1
                outer_annulus[x,y] = 1
            if sumsqu > radii and sumsqu < 2*radii:
                inner_annulus[x,y] = 1
                middle_annulus[x,y] = 0
                outer_annulus[x,y] = 1
            if sumsqu > 2*radii and sumsqu < 3*radii:
                inner_annulus[x,y] = 1
                middle_annulus[x,y] = 1
                outer_annulus[x,y] = 0
            else:
                inner_annulus[x,y] = 1
                middle_annulus[x,y] = 1
                outer_annulus[x,y] = 1

    return inner_annulus, middle_annulus, outer_annulus

def band_plot(image_r, image_g, image_u, cmap):
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

def hist(list1, list2):
    """
    Recieves two lists of seeing values, creates a bin histogram of each and
    displays tthe resulting plots.
    """
    figure, ax_array = plt.subplots(1,2)
    ax_array[0].hist(list1, range=(3,10))
    ax_array[1].hist(list2, range=(2,2.6))
    plt.show()