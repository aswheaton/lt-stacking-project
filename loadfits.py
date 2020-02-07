import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from scipy import optimize

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
                print(year, band, " matched ", filename)
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

def max_value_centroid(image_data, **kwargs):
    """
    Receives an image array and returns the coordinates of the brightest pixel
    in that image array.
    Args:
        image_data (2darray): The image array to be searched.
    Returns:
        (x_max,y_max) (tuple): pixel coordinates of the maximum value of the
            image array.
    """
    x_max, y_max = np.where(image_data == np.amax(image_data))
    # plt.imshow(cutout)
    # plt.scatter(x_max, y_max, s=2, c='red', marker='o')
    # plt.show()
    return((x_max[0], y_max[0]))

def annuli_mask(array, center, radii):
    """
    Receives and array and returns a tuple of three masked annuli from the
    input array about a given center.
    """
    return inner_annulus, middle_annulus, outer_annulus

def create_mask(image_data, **kwargs):
    # Offset image so that all values are positive
    image_data += np.abs(np.amin(image_data))
    mask = np.empty(image_data.shape)
    # Invalidate values based on the value of their neighbors (slow).
    if kwargs.get("condition") == "neighbors":
        for i in range(image_data.shape[0]):
            for j in range(image_data.shape[1]):
                try:
                    neighbors_sum = image_data[bc(i-1,j)]+image_data[bc(i+1,j)]+image_data[bc(i,j-1)]+image_data[bc(i,j+1)]
                    if image_data[i,j] >= neighbors_sum:
                        mask[i,j] = 0
                    else:
                        mask[i,j] = 1
                except IndexError:
                    mask[i,j] = 0
    # Invalidate values that fall below a certain threshold (fast).
    if kwargs.get("condition") == "threshold":
        max_value = np.amax(image_data)
        for i in range(image_data.shape[0]):
            for j in range(image_data.shape[1]):
                if image_data[i,j] <= 0.67 * max_value:
                    mask[i,j] = 0
                else:
                    mask[i,j] = 1
    return(mask)

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
    cutout += np.abs(np.amin(cutout))
    x_sum = np.sum(cutout, axis=0)
    y_sum = np.sum(cutout, axis=1)
    x_avg = np.average(range(x_sum.size), weights=x_sum)
    y_avg = np.average(range(y_sum.size), weights=y_sum)
    # plt.imshow(cutout)
    # plt.scatter(x_avg, y_avg, s=2, c='red', marker='o')
    # plt.show()
    if kwargs.get("floor") == True:
        return((int(np.floor(x_avg)), int(np.floor(y_avg))))
    else:
        return((x_avg, y_avg))

def hybrid_centroid(image_data, **kwargs):
    """
    Recieves an array of image data and returns the pixel coordinates of the
    centroid of the brightest star in the frame. Makes an initial guess at the
    position of the star by finding the maximum value in the array, then
    performs a weighted mean in two dimensions about the guess for finer accuracy.
    Args:
        image_data (2darray): array of image data containing reference star.
        size (int): the radius of the reference star, in pixels. Used to create
            cutout of appropriate size.
    Returns:
        (x_avg,y_avg) (tuple): pixel coordinates of the centroid of the
            brightest star in the image array.
    """
    proper_coords = degrees(("09:45:11.08","17:45:44.80"))
    # Attempt to invalidate pixels which may confuse the initial guess.
    if kwargs.get("mask") == True:
        masked_data = np.ma.array(image_data, mask=create_mask(image_data, condition="neighbors"))
        x_guess, y_guess = max_value_centroid(masked_data)
    elif kwargs.get("wcs") == True:
        x_guess, y_guess = wcs_offset(proper_coords, image_data)
    # Get the maximum value of the cutout as an initial guess.
    else:
        x_guess, y_guess = max_value_centroid(image_data)
    # Create a smaller cutout around the initial guess.
    cutout = np.array(image[0].data[x_guess-size:x_guess+size,y_guess-size:y_guess+size])
    x_avg, y_avg = weighted_mean_2D(cutout)
    # plt.imshow(cutout)
    # plt.scatter(x_avg, y_avg, s=2, c='red', marker='o')
    # plt.show()
    return((x_avg, y_avg))

def manual_centroid(image, **kwargs):
    """
    Allows the user to manually define the initial guess for centroiding by
    clicking on an imshow plot.

    Recieves an array of image 

def old_align(image_stack, **kwargs):
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

def align(images, **kwargs):
    """
    Recieves a list of image arrays containing a common object to use for
    alignment of the image stack. Returns a list of image arrays of different
    size, aligned, and with zero borders where the image has been shifted.
    Args:
        images (list): Frames to be aligned.
    Returns:
        aligned_images (list): new frames that have been aligned and can be
            stacked.
    """
    # Boolean, whether or not to mask images for hot pixels on the detector.
    mask = kwargs.get("mask")
    # Find the centroid of the reference star in each image.
    x_centroids, y_centroids = [], []
    for image in images:
        x_centroids.append(hybrid_centroid(image[0].data, size=50)[0])
        y_centroids.append(hybrid_centroid(image[0].data, size=50)[1])
    max_pos = (max(x_centroids), max(y_centroids))
    min_pos = (min(x_centroids), min(y_centroids))
    max_dif = (max_pos[0]-min_pos[0], max_pos[1]-min_pos[1])
    # Create new stack of aligned images using the centroid in each frame.
    aligned_images = []
    for image in images:
        aligned_image = np.zeros((image[0].data.shape[0]+max_dif[0], image[0].data.shape[1]+max_dif[1]))
        disp = (max_pos[0] - hybrid_centroid(image[0].data, size=50)[0], max_pos[1] - hybrid_centroid(image[0].data, size=50)[1])
        aligned_image[disp[0]:disp[0]+image[0].data.shape[0],disp[1]:disp[1]+image[0].data.shape[1]] = image[0].data
        aligned_images.append(aligned_image)
    return aligned_images

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
    proper_coords = degrees(("09:45:11.08","17:45:44.80"))

    for  year in range(2012,2018):
        for band in ["R", "G", "U"]:
            unaligned_images = load_fits(path="data/SDSSJ094511-P1-images/", year=str(year), band=band)
            aligned_images = align(unaligned_images, cutout=(x,y,dx,dy), centroid=hybrid_centroid)
            stacked_image = stack(aligned_images)
            rgu_images.append(stacked_image)

        plot(rgu_images[0], rgu_images[1], rgu_images[2], 'viridis')
        rgb(rgu_images[0][:1024,:1024], rgu_images[1][:1024,:1024], rgu_images[2][:1024,:1024])
        hist(seeing_pixels, seeing_arcsec)
main()
