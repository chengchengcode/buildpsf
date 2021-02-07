Here are two pieces of IDL codes to build psf from image. I use IRAC image as an example.

The main processs are:

1, make a rough catalog, select point source

2, check the point source, and reject the ones look bad

3, stack the point source after substract the background

I rebin the star image to align the brightest position in the center pixel.
I wrote these codes at 7 years before. I use them recently to build the PSF for the stacked IRAC image.

To check how good is the psf:

1, make growth curve and get the correction factors for different aperture photometry

2, do aperture photometry with different apertures

3, correct the aperture photometry by correction factor, and check if the values are consistent for point sources.

