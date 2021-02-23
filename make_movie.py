import os, sys
import numpy as np
import astropy.io.fits as pf
from astropy.coordinates import SkyCoord
hasistar = True
try:
    import istarmap   ### weird hack to make the progressbar work with starmap - see https://stackoverflow.com/questions/57354700/starmap-combined-with-tqdm/57364423#57364423
except:
    hasistar = False
    

from multiprocessing import Pool
import tqdm
import glob
  


class path:
    def __init__(self,ra,dec,imsize, framerate=24):
        self.framerate = int(framerate)  # frames per sec
        self.ra = np.array([ra])
        self.dec = np.array([dec])
        self.imsize = np.array([imsize])

    def panto(self,ra,dec,time):
        frames = int(time*self.framerate)
        
        ra0 = self.ra[-1]
        dec0 = self.dec[-1]
        imsize = self.imsize[-1]
        
        dsep = SkyCoord(ra0,dec0,unit='deg').separation(SkyCoord(ra,dec,unit='deg'))
        #th = SkyCoord(ra0,dec0,unit='deg').position_angle(SkyCoord(ra,dec,unit='deg')).value
        
        speed = dsep.to('arcsec').value / time
        print ('will pan {0:.1f} arcsec per second'.format(speed))
        ras = np.linspace(ra0,ra,frames)
        decs = np.linspace(dec0,dec,frames)
        self.add_path(ras, decs, np.ones_like(ras)*imsize)

    def goto(self,ra,dec,imsize=None):
        if imsize == None:
            imsize = self.imsize[-1]
        self.add_path(np.array([ra]),np.array([dec]),np.array([imsize]))
        
    def zoom(self,imsize1,time):
        frames = int(time*self.framerate)
        ra = self.ra[-1]
        dec = self.dec[-1]
        imsize0 = self.imsize[-1]
        imsizes = np.linspace(imsize0, imsize1, frames)
        
        self.add_path(ra*np.ones_like(imsizes), dec*np.ones_like(imsizes), imsizes)

    def dwell(self,time):
        '''
        time in sec
        '''
        frames = int(time*self.framerate)
        ra = self.ra[-1]
        dec = self.dec[-1]
        imsize = self.imsize[-1]
        self.add_path(ra*np.ones(frames), dec*np.ones(frames), imsize*np.ones(frames))
        
    def add_path(self,ra,dec,imsize):
        self.ra = np.hstack((self.ra, ra))
        self.dec = np.hstack((self.dec, dec))
        self.imsize = np.hstack((self.imsize, imsize))
        
        
def make_im(si,rai,deci,imsizei):
    '''
    call python script to make a single frame, number si at rai,deci with size imsizei
    '''
    cmd = 'python make_field_plot.py {im} frames/t{s:04d}.png --rms 1 --vmin 1.5 --vmax 150  --zoom_centre {ra:.5f},{dec:.5f} --zoom_size {size:.5f}'.format(im=fitsname, ra=rai,dec=deci,s=si,size=imsizei)
    #print(cmd)
    os.system(cmd)
    return cmd
        

if __name__ == '__main__':


    # define a path to travel through an image
    P = path(200.,55.,15.,framerate=24)
    P.dwell(1)
    P.zoom(4., 4.)  # zom to 4 deg fov in 3 s
    P.panto(210.7931525, 54.29116538, 30)
    P.zoom(2., 3.)
    P.dwell(1)
    P.goto(173.5785799, 49.07064301)
    P.dwell(1)
    P.zoom(4., 3.)
    P.panto(184.7855153, 47.27657956, 30)
    P.zoom(10., 2.)
    P.dwell(1)
    
    
    ra = P.ra
    dec = P.dec
    imsizes = P.imsize
    

    # make the frames
    Npoints = len(ra)
    print(('making {n} images'.format(n=Npoints)))
    print('with a framerate of {r:d} the movie will be {s:f} seconds'.format(r=P.framerate,s=Npoints/P.framerate))

    if os.path.isdir('frames'):
        print ('warning: output directory "frames" already exists')
        os.system('rm -rf frames/*')
    os.system('mkdir frames')

    fitsname = 'mosaic3.fits_snr.fits'

    
    ncpu = 12
    
    inputs = zip(range(Npoints),ra,dec,imsizes)
    with Pool(ncpu) as p:
        if hasistar:
            for _ in tqdm.tqdm(p.istarmap(make_im, inputs), total=Npoints):
                pass
        else: 
            p.starmap(make_im, inputs)

    # check all the frames are there (istarmap always seems to miss 1 or 2 of the first few?!?!?)
    tlist = sorted(glob.glob('frames/t*.png'))
    frames = np.array([int(t.replace('frames/t','').replace('.png','')) for t in tlist])
    missing =  np.array([fi  for fi in range(Npoints) if fi not in frames])
    print (missing)
    print (len(missing),'missing frames, trying to remake')
    
    inputs = zip(np.arange(Npoints)[missing],ra[missing],dec[missing],imsizes[missing])
    with Pool(ncpu) as p:
        p.starmap(make_im, inputs)
    
    # make the movie from the frames - note: need to have all the images (will complain if they're not consecutive)
    if os.path.exists('movie.mp4'):
        print ('warning, removing movie.mp4')
    os.system('rm movie.mp4')
    
    cmd = 'ffmpeg -f image2 -r {r:d} -i frames/t%04d.png movie.mp4'.format(r=P.framerate)
    print(cmd)
    os.system(cmd)

