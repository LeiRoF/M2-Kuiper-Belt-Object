
# File M2_CompuPhys_Data_Processing.py for Astronomical data processing course
#
# Example of Python code for measuring magnitudes in a fits image... Needs to be adapted !

def find_object(file, Xcent=1083, Ycent=660, r=8, r_in=15, r_out = 25, box_size = 80, plot = False, verbose = False):

    import matplotlib.pyplot as plt
    from astropy.io import fits

    hdulist = fits.open(file) # Define the image name to open
    if verbose: hdulist.info() # To get information on the image
    header = hdulist[0].header
    scidata = hdulist[0].data # transfer the image (pixels values) in scidata
    if verbose: print(scidata.shape)

    mjd = header['MJD-OBS']
    if verbose: print(mjd)

    exptime = header['EXPTIME']
    if verbose: print(exptime)

    #print(repr(header)) # to print all the descriptors in the header, if necessary
    #print(repr(header['MJD-OBS'])) # to print a given descriptor, if necessary
    #scidata /= exptime # divide the image by its exposure time, if necessary

    # r=8 # radius for measuring target flux
    # r_in=15 # inner radius for the annulus used for measuring sky background
    # r_out=25 # outer radius for the annulus used for measuring sky background

    import numpy as np
    from photutils import CircularAperture, CircularAnnulus, aperture_photometry
    from astropy.visualization import simple_norm

    # Xcent=1083 # approximate target coordinates (x)
    # Ycent=660 # approximate target coordinates (y)

    from photutils.centroids import centroid_com, centroid_1dg, centroid_2dg, centroid_quadratic

    D = box_size//2

    data = scidata[Ycent-D:Ycent+D,Xcent-D:Xcent+D] # subimage centered on the target
    x1, y1 = centroid_quadratic(data) # for computing accurately the position of maximum of brightness
    if verbose: print('resultats centroid =%.3f'%(x1),',%.3f'%(y1))

    for _ in range(1): # avoid problem with syntaxic coloration caussed by simple_norm()
        norm = simple_norm(data, 'sqrt', percent=99)

    if plot:
        plt.imshow(data, norm=norm, origin='lower')
        positions = [(x1, y1)]
        aperture = CircularAperture(positions, r)
        annulus_aperture = CircularAnnulus(positions, r_in, r_out)
        aperture.plot(color='white', lw=2)
        annulus_aperture.plot(color='white', lw=2)
        phot_table = aperture_photometry(data, aperture)
        phot_table['aperture_sum'].info.format = '%.8g' # for consistent table output
        Sum_target_raw=(phot_table['aperture_sum'])
        phot_table = aperture_photometry(data, annulus_aperture)
        phot_table['aperture_sum'].info.format = '%.8g' # for consistent table output
        Sky_background=(phot_table['aperture_sum'])
        from math import log10
        Sum_target=Sum_target_raw-Sky_background/(r_out*r_out-r_in*r_in)*(r*r)
        Magnitude=-2.5*log10(Sum_target)
    
    if verbose:
        print(phot_table)
        print('\nSum_target_raw=%.3f\n'%(Sum_target_raw))
        print(phot_table)
        print('\nSky_background=%.3f\n'%(Sky_background))
        print('Target flux without sky=%.3f\n'%(Sum_target))
        print('instrumental magnitude=%.3f\n'%(Magnitude))

    if plot:
        plt.title(file)
        plt.show()

    hdulist.close()
    return [x1+Xcent-box_size//2, y1+Ycent-box_size//2]