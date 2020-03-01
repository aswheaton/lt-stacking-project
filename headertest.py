"""
CCDXIMSI=                 1024 / [pixels] Imaging pixels
CCDYIMSI=                 1024 / [pixels] Imaging pixels
CCDXBIN =                    2 / [pixels] X binning factor
CCDYBIN =                    2 / [pixels] Y binning factor
CCDXPIXE=            0.0000135 / [m] Size of pixels, in X:13.5um
CCDYPIXE=            0.0000135 / [m] Size of pixels, in Y:13.5um
CCDSCALE=                0.279 / [arcsec/binned pixel] Scale of binned image on
RA      = '09:45:10.866'       / World coordinate at the reference pixel
DEC     = '+17:45:48.17'       / World coordinate at the reference pixel
RADECSYS= 'FK5     '           / [FK4, FK5] Fundamental coord system
L1SEEING=             4.350539 / [pixels] frame seeing in pixels
L1SEESEC=             1.213801 / [arcsec] frame seeing in arcsec
SEEING  =             4.350539 / [pixels] frame seeing in pixels
CTYPE1  = 'RA---TAN'           / WCS projection
CTYPE2  = 'DEC--TAN'           / WCS projection
CRPIX1  =                 512. / WCS reference pixel coordinate
CRPIX2  =                 512. / WCS reference pixel coordinate
CRVAL1  =        146.295274426 / [degrees] World coordinate at the ref pix
CRVAL2  =          17.76338016 / [degrees] World coordinate at the ref pix
WRA     = '9:45:11.075'        / Original RA value from TCS before WCS fit
WDEC    = '+17:45:44.89'       / Original DEC value from TCS before WCS fit
WRADECSY= 'FK5     '           / [FK4,FK5] Original Fundamental system before WC
WROTSKY =             270.0001 / [deg] Original sky PA from TCS before WCS fit
EPOCH   =                2000. / Epoch of coordinate
CDELT1  =       -7.7508973E-05 / [degrees/pixel]
CDELT2  =        7.7508973E-05 / [degrees/pixel]
CROTA1  =            90.361783 / [degrees]
CROTA2  =            90.361783 / [degrees]
"""

from astropy.io import fits

image = fits.open("data/fits/20120312_38_G100.fits")

image.info()

for key in image[0].header:
    print(key, image[0].header[key])
