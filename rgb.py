"""
https://docs.astropy.org/en/stable/visualization/rgb.html
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import make_lupton_rgb
image_2012_r = fits.open("data/fits/20120312_39_R100.fits")[0].data / 1.5
image_2012_g = fits.open("data/fits/20120312_38_G100.fits")[0].data
image_2012_b = fits.open("data/fits/20120312_40_U300.fits")[0].data
image_2017_r = fits.open("data/fits/20170314_20_R.fits")[0].data / 3.0
image_2017_g = fits.open("data/fits/20170314_19_G.fits")[0].data / 3.0
image_2017_b = fits.open("data/fits/20170314_21_U.fits")[0].data
image_2012 = make_lupton_rgb(image_2012_r, image_2012_g, image_2012_b, Q=10, stretch=1000.)
plt.imshow(image_2012)
plt.show()
image_2017 = make_lupton_rgb(image_2017_r, image_2017_g, image_2017_b, Q=10, stretch=5000.)
plt.imshow(image_2017)
plt.show()
