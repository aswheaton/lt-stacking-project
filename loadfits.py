import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np

hdul = fits.open("data/fits/20120312_38_G100.fits")
hdul.info()
plt.imshow(hdul[0].data)
plt.show()
