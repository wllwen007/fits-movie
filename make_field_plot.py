# W Williams
# utility for making movie frame images of fits files
#

import matplotlib as mpl
mpl.use('Agg')
import aplpy as ap
import pyregion
from astropy.visualization import SqrtStretch,ImageNormalize

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import os
import argparse

import warnings
warnings.filterwarnings("ignore", category=FutureWarning) 


def flatten(f,ra,dec,x,y,size,size2=None,hduid=0,channel=0,freqaxis=3,verbose=True):
    ''' 
    Flatten a fits file so that it becomes a 2D image. Return new header and
    data
    This version also makes a sub-image of specified size.
    modified to require full image in FoV
    
    from: https://github.com/mhardcastle/lotss-catalogue/tree/master/utils/subim.py
    '''

    naxis=f[hduid].header['NAXIS']
    if naxis<2:
        raise RuntimeError('Can\'t make map from this')

    if size2 is None:
        size2 = size

    if verbose:
        print(f[hduid].data.shape)
    ds=f[hduid].data.shape[-2:]
    by,bx=ds
    if (x<0) or (x>bx):
        raise RuntimeError('coordinates not in image')
    if (y<0) or (y>by):
        raise RuntimeError('coordinates not in image')
        
    xmin=int(x-size)
    if xmin<0:
        raise RuntimeError('zoom too large for image')
    xmax=int(x+size)
    if xmax>bx:
        raise RuntimeError('zoom too large for image')
    ymin=int(y-size2)
    if ymin<0:
        raise RuntimeError('zoom too large for image')
    ymax=int(y+size2)
    if ymax>by:
        raise RuntimeError('zoom too large for image')
    
    if ymax<=ymin or xmax<=xmin:
        # this can only happen if the required position is not on the map
        print(xmin,xmax,ymin,ymax)
        raise RuntimeError('Failed to make subimage!')

    w = WCS(f[hduid].header)
    wn=WCS(naxis=2)
    
    wn.wcs.crpix[0]=w.wcs.crpix[0]-xmin
    wn.wcs.crpix[1]=w.wcs.crpix[1]-ymin
    wn.wcs.cdelt=w.wcs.cdelt[0:2]
    try:
        wn.wcs.pc=w.wcs.pc[0:2,0:2]
    except AttributeError:
        pass # pc is not present
    wn.wcs.crval=w.wcs.crval[0:2]
    wn.wcs.ctype[0]=w.wcs.ctype[0]
    wn.wcs.ctype[1]=w.wcs.ctype[1]
    
    header = wn.to_header()
    header["NAXIS"]=2

    slice=[]
    for i in range(naxis,0,-1):
        if i==1:
            slice.append(np.s_[xmin:xmax])
        elif i==2:
            slice.append(np.s_[ymin:ymax])
        elif i==freqaxis:
            slice.append(channel)
        else:
            slice.append(0)
    if verbose:
        print(slice)

    header['HISTORY'] = 'flatten'
    hdu=fits.PrimaryHDU(f[hduid].data[slice],header)
    copy=('EQUINOX','EPOCH','BMAJ','BMIN','BPA')
    for k in copy:
        r=f[hduid].header.get(k)
        if r is not None:
            hdu.header[k]=r
    if 'TAN' in hdu.header['CTYPE1']:
        hdu.header['LATPOLE']=f[hduid].header['CRVAL2']
    hdulist=fits.HDUList([hdu])
    return hdulist


def extract_subim(filename,ra,dec,size,size2=None,hduid=0,verbose=False):
    if verbose:
        print('Opening',filename)
    if size2 is None:
        size2 = size
    orighdu=fits.open(filename)
    psize=int(size/abs(orighdu[hduid].header['CDELT1']))
    psize2=int(size2/abs(orighdu[hduid].header['CDELT2']))
    if verbose:
        print('pix size is',psize,size)
        print('pix size2 is',psize2,size2)
    ndims=orighdu[hduid].header['NAXIS']
    pvect=np.zeros((1,ndims))
    lwcs=WCS(orighdu[hduid].header)
    pvect[0][0]=ra
    pvect[0][1]=dec
    imc=lwcs.wcs_world2pix(pvect,0)
    x=imc[0][0]
    y=imc[0][1]
    hdu=flatten(orighdu,ra,dec,x,y,psize,psize2,hduid=hduid,verbose=verbose)
    return hdu


def plot_image(fitsname, imagename, rms=np.nan, vmin=-3, vmax=25, ra_zoom=None, dec_zoom=None, d_zoom=None, AR=16/9.):
    '''
    make an image called imagename from cutout from fits file fitsname at coords ra_zoom, dec_zoom and dimensions d_zoom x d_zoom/AR
    '''
    
    hdu = extract_subim(fitsname, ra_zoom, dec_zoom, d_zoom, d_zoom/AR)

    dat = hdu[0].data
    figsize= (10,10./AR)
    f1 = plt.figure(figsize = figsize)
    
    norm = ImageNormalize(vmin=vmin*rms, vmax=vmax*rms, stretch=SqrtStretch())    
    plt.imshow(dat, origin='lower', cmap=plt.cm.cubehelix, norm=norm)

    
    plt.subplots_adjust(left=0.0, bottom=0.0, top=1.0, right=1.0)
    
    plt.savefig(imagename)


    return

def main():
    
    parser = argparse.ArgumentParser() #prog='plot_solutions_all_stations.py',usage='[options] <parmdb> <imageroot> ')
    parser.add_argument('--rms', default=None, help="RMS to use")
    parser.add_argument('--vmin', default=-3., help="Min scaling values (in sigma units)")
    parser.add_argument('--vmax', default=25., help="Max scaling values (in sigma units)")
    parser.add_argument('--zoom_centre', default=None, help="Centre of cutout (ra,dec in degrees)")
    parser.add_argument('--zoom_size', default=1., help="size of cutout (in degrees) (default: 1 deg)")
    parser.add_argument('fits_image', help="Name of fits image")
    parser.add_argument('output_image', help="Name of output image")
    
    args = parser.parse_args()
    
    fits_image = args.fits_image
    output_image = args.output_image
    rms = args.rms
    if rms is not None:
        rms = float(rms)
    else:
        rms = np.nan
    zoom_centre = args.zoom_centre
    if zoom_centre is not None:
        c = zoom_centre.split(',')
        ra_zoom = float(c[0])
        dec_zoom  = float(c[1])
    else:
        ra_zoom = None
        dec_zoom = None
    d_zoom = float(args.zoom_size)
    vmin = float(args.vmin)
    vmax = float(args.vmax)
    
    plot_image(fits_image, output_image, rms=rms, vmin=vmin, vmax=vmax, ra_zoom=ra_zoom, dec_zoom=dec_zoom, d_zoom=d_zoom)



if __name__ == "__main__":
    main()
