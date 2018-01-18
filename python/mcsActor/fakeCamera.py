import numpy
import time
import pyfits

class Camera(object):
    pass

class FakeCamera(Camera):
    def __init__(self):
        # Should read from config file....
        self.imageSize = (4096, 4096)
        self.biasLevel = 100
        self.readNoise = 5.0
        self.name = 'numpy_fake'
        
    def _readoutTime(self):
        return 0.5

    def _exposureOverheadTime(self):
        return 0.1

    def sendStatusKeys(self, cmd):
        """ Send our status keys to the given command. """ 

        cmd.inform('cameraName=%s; readNoise=%0.2f' % (self.name, self.readNoise))

    def expose_fake(self, cmd, filename, get_arc):

        """

        Read an image file (and optionally an arc file) from disk, to mimic
        a real exposure.

        Args: 
            cmd     - a Command object to report to. Ignored if None.
            filename - FITS filename (minus FITS extension) of the file to read into
            get_arc - flag for arc image. If get_arc=1, looks for filename_arc.fits

        Returns 
            the image
            the arc image (optional)

        """
        
        image=pyfits.getdata(filename+".fits",1)

        if(get_arc==1):
            arc_image=pyfits.getdata(filename+"_arc.fits",1)
            return image,arc_iamge
        else:
            return image
     
    def expose(self, cmd, expTime, expType):
        """ Generate an 'exposure' image. We don't have an actual camera, so generate some 
        plausible image. 

        Args:
           cmd     - a Command object to report to. Ignored if None.
           expTime - the fake exposure time. 
           expType - ("bias", "dark", "object", "flat", "test")
           
        Returns:
           - the image.

        Keys:
           exposureState
        """

        if not expType:
            expType = 'test'
        if cmd:
            cmd.inform('exposureState="exposing"')
        if expType not in ('bias', 'test') and expTime > 0:
            time.sleep(expTime + self._exposureOverheadTime())

        if cmd:
            cmd.inform('exposureState="reading"')
        
        f = pyfits.open('/home/chyan/mhs/data/mcs/schmidt_fiber_snr400_rmod71.fits')      
        image = f[0].data
        #image = numpy.random.normal(self.biasLevel, 
        #                            scale=self.readNoise, 
        #                            size=self.imageSize).astype('u2')

        if expType != 'test':
            time.sleep(self._readoutTime())
        return image
        
        
