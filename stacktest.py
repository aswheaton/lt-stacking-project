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
x_res, y_res = np.linspace(0, 20, 100), np.linspace(0, 20, 100)
x_res, y_res  = np.meshgrid(x_res, y_res)
stack = np.zeros((100,100))

for image in unaligned_images:
    cutout = wcs_cutout(image, size=10, proper_coords=proper_coords)
    initial_guess = get_gauss_guess(cutout)
    params, covariance = gaussian_fit(data=cutout, guess=initial_guess)
    fitted_data = gaussian_2D((x_res, y_res), *params).reshape(100, 100)
    plt.imshow(cutout, cmap='viridis', origin='lower', norm=LogNorm())
    plt.contour(x_res, y_res, fitted_data, 7, colors='r')
    plt.show()
    # Re-center the gaussian fit on the center of the cutout.
    c_params = (params[0], 10.0, 10.0, params[3], params[4], params[5], params[6])
    stack += gaussian_2D((x_res,y_res), *c_params).reshape(100,100)
plt.imshow(stack, cmap='viridis', origin='lower', norm=LogNorm())
plt.show()
plt.clf()
# plot(rgu_images[0], rgu_images[1], rgu_images[2], 'viridis')
# rgb(rgu_images[0][:1024,:1024], rgu_images[1][:1024,:1024], rgu_images[2][:1024,:1024])
# hist(seeing_pixels, seeing_arcsec)
