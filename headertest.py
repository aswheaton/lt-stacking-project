from astropy.io import fits

image = fits.open("data/fits/20120312_38_G100.fits")

image.info()

for key in image[0].header:
    print(key, image[0].header[key])
