from astropy.io import fits
import numpy as np

data = fits.open("data/fits/20120312_38_G100.fits")
data.info()
print(np.array(data[0]))
