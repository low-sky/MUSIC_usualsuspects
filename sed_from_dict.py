import numpy as np
import pylab as pl

band_waves = {0: 1.98, 1: 1.3, 2: 1.0, 3: 0.86}

def sed_from_dict(data, TestRadius=1.1, InnerApp=20, OuterApp=23):
    """
    TestRadius : float
        Pixel search radius to examine for peak location...
    """
    shape = data[0].shape
    cube = np.zeros((len(data),shape[0],shape[1]))
    for band in band_waves:
        cube[band,:,:] = data[band]

    yy,xx = np.indices(shape)

    # original radius, centered at center...
    xcen,ycen=np.asarray(shape)/2
    rr = ((yy-ycen)**2 + (xx-xcen)**2)**0.5

    # recenter on peak
    xcen,ycen = np.unravel_index((np.argmax((rr<TestRadius)*cube[1,:,:])),shape)
    rr = ((yy-ycen)**2 + (xx-xcen)**2)**0.5

    rr1,bandidx = np.meshgrid(rr,xrange(4))
    rr1.shape = cube.shape
    bandidx.shape = cube.shape

    # not used counts = np.nansum(np.nansum((rr1 < TestRadius)*(cube>0),axis=1),axis=1)
    fluxes = np.nansum(np.nansum(cube*(rr1 < TestRadius)*(cube>0),axis=1),axis=1)

    medians = np.zeros(4)
    mads = np.zeros(4)

    for i in xrange(4):
        indices = (rr1>InnerApp)*(rr1<OuterApp)*(cube==cube)*(bandidx == i)
        medians[i] = np.median(cube[indices])
        mads[i] = np.median(np.abs(cube[indices]-medians[i]))/0.6745


    return fluxes,medians,mads

def plot_sed(fluxes, backgrounds, errors, **kwargs):
    pl.errorbar(band_waves.values(),fluxes-backgrounds,yerr=errors,marker='s', **kwargs)
    pl.xlabel('$\lambda$ (mm)')
    pl.ylabel('mJy/beam')
