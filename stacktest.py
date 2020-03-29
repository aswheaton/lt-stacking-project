import sys
from fits_utils import *

# Empty list for collecting stacked images in three bands.
rgu_images = []
# Define the cutout region containing the reference object.
x, y, dx, dy = 450, 450, 75, 60
# Define the proper location of the object, for alignment.
proper_coords = degrees(("09:45:11.08","17:45:44.80"))

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

bins=[0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75]

print("Max: {}".format(max(seeing_vals)))
print("Mean: {}".format(np.mean(seeing_vals)))
print("StDev: {}".format(np.std(seeing_vals)))
plt.hist(seeing_vals, bins=bins)
plt.axvline(np.mean(seeing_vals)-np.std(seeing_vals), color='r', linestyle='dashed', linewidth=1)
plt.axvline(np.mean(seeing_vals), color='b', linestyle='dashed', linewidth=1)
plt.axvline(np.mean(seeing_vals)+np.std(seeing_vals), color='r', linestyle='dashed', linewidth=1)
plt.axvline(np.mean(seeing_vals)+2*np.std(seeing_vals), color='r', linestyle='dashed', linewidth=1)
plt.axvline(np.mean(seeing_vals)+3*np.std(seeing_vals), color='r', linestyle='dashed', linewidth=1)
plt.title("Astronomical Seeing in {}-band".format(band))
plt.xlim([min(bins),max(bins)])
plt.xticks(bins)
plt.savefig("report/img/seeing_hist_{}_band.eps".format(band),bbox_inches="tight", pad_inches=0)
plt.show()
plt.clf()

# for image in unaligned_images:
#     obj_ra, obj_dec = wcs_centroid(proper_coords, image)
#     plt.imshow(image["data"][475:575,450:550], cmap='viridis', origin='lower', norm=LogNorm())
#     plt.scatter(obj_dec-475.0, obj_ra-475.0, s=2, c='red', marker='o')
#     plt.show()
#     plt.clf()

# aligned_images = align(unaligned_images, centroid=manual_centroid)
# stacked_image = stack(aligned_images, correct_exposure=False)
# plt.imshow(stacked_image["data"][475:575,450:550], cmap='viridis', origin='lower', norm=LogNorm())

# plt.savefig("plots/manual_centroid_{}_stack.eps".format(band),
#             bbox_inches="tight", pad_inches="0", dpi=1200, format="eps"
#             )
# plt.savefig("plots/manual_centroid_{}_stack.png".format(band),
#             bbox_inches="tight", pad_inches=0, dpi=1200, format="png"
#             )

# plt.show()

# plot(rgu_images[0], rgu_images[1], rgu_images[2], 'viridis')
# rgb(rgu_images[0][:1024,:1024], rgu_images[1][:1024,:1024], rgu_images[2][:1024,:1024])
# hist(seeing_pixels, seeing_arcsec)
