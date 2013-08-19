from astropy.io import fits
import idlsave
import itertools
import scipy.ndimage
from astropy.nddata.convolution import make_kernel,convolve
import numpy as np

data_w49_1   = {k:idlsave.read('w49/130819_ob5_band%ii_clean_music_20130815_map.sav' % k) for k in xrange(4)}
data_w49_2   = {k:idlsave.read('w49/130819_ob6_band%ii_clean_music_20130815_map.sav' % k) for k in xrange(4)}
data_w51_1   = {k:idlsave.read('w51_bgps/130819_ob1_band%ii_clean_music_20130815_map.sav' % k) for k in xrange(4)}
data_w51_2   = {k:idlsave.read('w51_bgps/130819_ob2_band%ii_clean_music_20130815_map.sav' % k) for k in xrange(4)}
data_sgrb2_1 = {k:idlsave.read('sgr_b2/130819_ob3_band%ii_clean_music_20130815_map.sav' % k) for k in xrange(4)}
data_sgrb2_2 = {k:idlsave.read('sgr_b2/130819_ob4_band%ii_clean_music_20130815_map.sav' % k) for k in xrange(4)}

datasets = dict(zip(['W49_1','W49_2','W51_1','W51_2','SgrB2_1','SgrB2_2'],(data_w49_1, data_w49_2, data_w51_1, data_w51_2, data_sgrb2_1, data_sgrb2_2)))

for jj,(k,data) in enumerate(datasets.iteritems()):

    band3 = data[3].mapstruct.map[0]
    yy1,xx1 = grid1 = np.indices(band3.shape)

    pixscale = data[3].mapstruct['OMEGA_PIX_AM']**0.5 / 60.

    header = fits.Header()
    header['CDELT1'] = pixscale
    header['CDELT2'] = pixscale
    ffile = fits.PrimaryHDU(data=band3, header=header)

    for ii in xrange(4):
        m = data[ii].mapstruct.map[0]
        ratio = m.shape[0]/float(band3.shape[0])
        newm = scipy.ndimage.map_coordinates(np.nan_to_num(m), grid1*ratio)
        bads = scipy.ndimage.map_coordinates(np.array(m!=m,dtype='float'), grid1*ratio)
        newm[bads>0.5] = np.nan

        beamsize_delta = (np.abs(data[ii].mapstruct['OMEGA_BEAM_AM']-data[0].mapstruct['OMEGA_BEAM_AM'])/np.pi)**0.5
        am_per_pix = data[0].mapstruct['OMEGA_PIX_AM']**0.5
        kernelwidth = beamsize_delta/am_per_pix
        if kernelwidth > 0:
            kernel = make_kernel.make_kernel(band3.shape, kernelwidth=kernelwidth)
            newm = convolve.convolve_fft(newm, kernel, interpolate_nan=True)

        newm /= data[3].mapstruct['OMEGA_BEAM_AM']/data[ii].mapstruct['OMEGA_BEAM_AM']

        kernel = make_kernel.make_kernel(band3.shape, kernelwidth=5)
        smm = convolve.convolve_fft(newm,kernel, interpolate_nan=True)

        ffile.data = newm

        ffile.writeto("%s_Band%i_smooth.fits" % (k,ii))

        ffile.data = newm - smm

        ffile.writeto("%s_Band%i_smooth_unsharp.fits" % (k,ii))
