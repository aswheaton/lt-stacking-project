import sys
from fits_utils import *

# Empty list for collecting stacked images in three bands.
rgu_images = []
# Define the cutout region containing the reference object.
x, y, dx, dy = 450, 450, 75, 60
# Define the proper location of the object, for alignment.
proper_coords = degrees(("09:45:11.08","17:45:44.78"))
print(proper_coords)

band = str(sys.argv[1])
unaligned_images = []

for year in range(2012,2019):
    unaligned_images += load_fits(path="data/fits/", year=str(year), band=band)

for image in unaligned_images:
    # Get the rotation of the image in the stack.
    quarter_rotations_1 = int(np.around(image["rot1"] / 90.0))
    quarter_rotations_2 = int(np.around(image["rot2"] / 90.0))
    if quarter_rotations_1 == quarter_rotations_2:
        q_rot = quarter_rotations_1
    else:
        print("Error: Image rotations do not match!")
        break
    # Correct rotation, if needed.
    image["data"] = np.rot90(image["data"], q_rot)

for image in unaligned_images:
    size_1, size_2 = image["data"].shape
    if size_1 != 1024 or size_2 != 1024:
        image["data"] = image["data"][512:1536,512:1536]
    # plt.imshow(image["data"], cmap='viridis', origin='lower', norm=LogNorm())
    # plt.show()
    # plt.clf()

seeing_vals = []
for image in unaligned_images:
    if image["seeing"] > 3.0:
        pass
    else:
        seeing_vals.append(image["seeing"])
seeing_filtered_images = []
for image in unaligned_images:
    if image["seeing"] <= np.mean(seeing_vals) + np.std(seeing_vals):
        seeing_filtered_images.append(image)
unaligned_images = seeing_filtered_images
# print("Max: {}".format(max(seeing_vals)))
# print("Mean: {}".format(np.mean(seeing_vals)))
# print("StDev: {}".format(np.std(seeing_vals)))
# bins=[0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75]
# plt.hist(seeing_vals, bins=bins)
# plt.axvline(np.mean(seeing_vals)-np.std(seeing_vals), color='r', linestyle='dashed', linewidth=1)
# plt.axvline(np.mean(seeing_vals), color='b', linestyle='dashed', linewidth=1)
# plt.axvline(np.mean(seeing_vals)+np.std(seeing_vals), color='r', linestyle='dashed', linewidth=1)
# plt.axvline(np.mean(seeing_vals)+2*np.std(seeing_vals), color='r', linestyle='dashed', linewidth=1)
# plt.axvline(np.mean(seeing_vals)+3*np.std(seeing_vals), color='r', linestyle='dashed', linewidth=1)
# plt.title("Astronomical Seeing in {}-band".format(band))
# plt.xlim([min(bins),max(bins)])
# plt.xticks(bins)
# plt.savefig("report/img/seeing_hist_{}_band.eps".format(band),bbox_inches="tight", pad_inches=0)
# plt.show()
# plt.clf()

# center, bb_low, bb_high = 512, 350, 650
# for image in unaligned_images:
#     obj_ra, obj_dec = wcs_centroid(image, proper_coords=proper_coords, correction_factor=(0,0))
#     plt.imshow(image["data"][bb_low:bb_high,bb_low:bb_high], cmap='viridis', origin='lower', norm=LogNorm())
#     plt.scatter(center-bb_low, center-bb_low, s=8, c='red', marker='o')
#     plt.scatter(obj_ra-bb_low, obj_dec-bb_low, s=2, c='red', marker='o')
#     plt.show()
#     plt.clf()

aligned_images = align(unaligned_images, centroid=wcs_centroid, proper_coords=proper_coords)
stacked_image = stack(aligned_images, correct_exposure=False)
plt.imshow(stacked_image["data"][553:573,565:585], cmap='viridis', origin='lower', norm=LogNorm())
plt.title("J094511, {}-band".format(band))
# plt.savefig("report/img/wcs_centroid_{}_stack.eps".format(band),bbox_inches="tight", pad_inches=0)
plt.show()

cropped_stack = np.array(stacked_image["data"][553:573,565:585])
initial_guess = get_gauss_guess(cropped_stack)
params, covariance = gaussian_fit(data=cropped_stack, guess=initial_guess)

# Plot a contour map over the stacked object image.
x = np.linspace(0, cropped_stack.shape[0], cropped_stack.shape[0])
y = np.linspace(0, cropped_stack.shape[1], cropped_stack.shape[1])
x, y  = np.meshgrid(x, y)
fitted_data = gaussian_2D((x, y), *params).reshape(20, 20)
plt.imshow(cropped_stack, cmap='viridis', origin='lower', norm=LogNorm())
plt.contour(x, y, fitted_data, 7, colors='r')
plt.title("Gaussian Fit on J094511, {}-band".format(band))
# plt.savefig("report/img/gauss_fit_wcs_{}_stack.eps".format(band),bbox_inches="tight", pad_inches=0)
plt.show()
plt.clf()

"""
Now fit a gaussian to the object in each unstacked frame, resample, and stack.
"""
x, y = np.linspace(0, 20, 20), np.linspace(0, 20, 20)
x, y = np.meshgrid(x, y)
x_res, y_res = np.linspace(0, 20, 1000), np.linspace(0, 20, 1000)
x_res, y_res  = np.meshgrid(x_res, y_res)
stack = np.zeros((1000,1000))

for image in unaligned_images:
    cutout = wcs_cutout(image, size=10, proper_coords=proper_coords)
    initial_guess = get_gauss_guess(cutout)
    params, covariance = gaussian_fit(data=cutout, guess=initial_guess)
    fitted_data = gaussian_2D((x_res, y_res), *params).reshape(1000, 1000)
    # plt.imshow(cutout, cmap='viridis', origin='lower', norm=LogNorm())
    # plt.contour(x_res, y_res, fitted_data, 7, colors='r')
    # plt.show()
    # Re-center the gaussian fit on the center of the cutout.
    c_params = (params[0], 10.0, 10.0, params[3], params[4], params[5], params[6])
    stack += gaussian_2D((x_res,y_res), *c_params).reshape(1000,1000)

# Crude background subtraction.
stack -= np.min(stack)
# Fit a final gaussian to the stacked gaussians and plot as a contour map.
initial_guess = get_gauss_guess(stack)
params, covariance = gaussian_fit(data=stack, guess=initial_guess)
fitted_data = gaussian_2D((x_res, y_res), *params).reshape(1000, 1000)
print("Amp: {}\nx0,y0: {}, {}\nsigma_x: {}\n sigma_y: {}\ntheta: {}\noffset: {}".format(*params))
plt.imshow(stack, cmap='viridis', origin='lower', norm=LogNorm())
x_res, y_res = np.linspace(0, 1000, 1000), np.linspace(0, 1000, 1000)
plt.contour(x_res, y_res, stack, 7, colors='r')
plt.show()
plt.clf()

steps = 49
radii = np.linspace(50.0, 500.0, steps)
total_counts = np.zeros(steps)
counts_per_area = np.zeros(steps)
annulus_area = np.zeros(steps)
annulus_counts = np.zeros(steps)

for i in range(steps):
    counts = 0.0
    for row in range(1000):
        for col in range(1000):
            p_radius = np.sqrt((row-500)**2 + (col-500)**2)
            if p_radius < radii[i]:
                counts += stack[row,col]
    total_counts[i] = counts
    annulus_counts[i] = counts - total_counts[i-1]
    annulus_area[i] = np.pi*(radii[i]**2 - (radii[i]-50.0)**2)
    counts_per_area[i] = annulus_counts[i] / annulus_area[i]

plt.plot(radii, total_counts, 'r-')
plt.title("Counts in Aperture")
plt.savefig("report/img/J094511_{}_aperture_sum.eps".format(band), bbox_inches="tight", pad_inches=0)
plt.savefig("report/img/J094511_{}_aperture_sum.png".format(band), bbox_inches="tight", pad_inches=0)
plt.show()
plt.plot(radii, annulus_counts, 'g-')
plt.title("Counts in Annulus")
plt.savefig("report/img/J094511_{}_annulus_sum.eps".format(band), bbox_inches="tight", pad_inches=0)
plt.savefig("report/img/J094511_{}_annulus_sum.png".format(band), bbox_inches="tight", pad_inches=0)
plt.show()
plt.plot(radii, annulus_area, 'b-')
plt.title("Area of Annulus")
plt.savefig("report/img/J094511_{}_annulus_area.eps".format(band), bbox_inches="tight", pad_inches=0)
plt.savefig("report/img/J094511_{}_annulus_area.png".format(band), bbox_inches="tight", pad_inches=0)
plt.show()
plt.plot(radii, counts_per_area, 'r-')
plt.title("Flux in Annulus (J094511 Profile)")
plt.savefig("report/img/J094511_{}_annulus_flux.eps".format(band), bbox_inches="tight", pad_inches=0)
plt.savefig("report/img/J094511_{}_annulus_flux.png".format(band), bbox_inches="tight", pad_inches=0)
plt.show()
# plot(rgu_images[0], rgu_images[1], rgu_images[2], 'viridis')
# rgb(rgu_images[0][:1024,:1024], rgu_images[1][:1024,:1024], rgu_images[2][:1024,:1024])
# hist(seeing_pixels, seeing_arcsec)
