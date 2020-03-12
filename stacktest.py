from fits_utils import *

# Empty list for collecting stacked images in three bands.
rgu_images = []
# Define the cutout region containing the reference object.
x, y, dx, dy = 450, 450, 75, 60
# Define the proper location of the object, for alignment.
proper_coords = degrees(("09:45:11.08","17:45:44.80"))

band = "G"
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
    plt.imshow(image["data"], cmap='viridis', origin='lower', norm=LogNorm())
    plt.show()
    plt.clf()


# aligned_images = align(unaligned_images, centroid=manual_centroid)
# stacked_image = stack(aligned_images, correct_exposure=False)
# plt.imshow(stacked_image["data"][475:575,450:550], cmap='viridis', origin='lower', norm=LogNorm())
# plt.show()

plot(rgu_images[0], rgu_images[1], rgu_images[2], 'viridis')
rgb(rgu_images[0][:1024,:1024], rgu_images[1][:1024,:1024], rgu_images[2][:1024,:1024])
hist(seeing_pixels, seeing_arcsec)
