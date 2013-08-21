import pylab as pl
import idlsave

def load_data(obsname):
    """ load data file given an obs date/number """
    filename = obsname+"_band%ii_clean_music_20130815_map.sav"
    data = {k:idlsave.read(filename % k, verbose=False) for k in xrange(4)}
    return data

def viewer(data, cb=False, **kwargs):
    for ii,(k,v) in enumerate(data.iteritems()):
        pl.subplot(2,2,ii+1)
        pl.imshow(v.mapstruct.map[0], **kwargs)
        if cb:
            pl.colorbar()

def dictviewer(data, **kwargs):
    for ii,(k,v) in enumerate(data.iteritems()):
        pl.subplot(2,2,ii+1)
        pl.imshow(v, **kwargs)
