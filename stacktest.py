from fits_utils import *

band = "G"
unaligned_images = []

for year in range(2012,2019):
    unaligned_images += load_fits(path="data/SDSSJ094511-P1-images/",
                                  year=str(year), band=band
                                  )

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
