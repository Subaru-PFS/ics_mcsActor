from builtins import object
import numpy
import time
import subprocess as sub
import astropy.io.fits as pyfits


class Camera(object):
    pass

class mcsCamera(Camera):
    def __init__(self):
        # Should read from config file....
        self.imageSize = (8960, 5778)
        self.biasLevel = 100
        self.readNoise = 5.0
        self.name = 'Canon_50M'
        self.expTime = 0
        
    def _readoutTime(self):
        return 0.5

    def _exposureOverheadTime(self):
        return 0.1

    def sendStatusKeys(self, cmd):
        """ Send our status keys to the given command. """ 

        cmd.inform('cameraName=%s; readNoise=%0.2f' % (self.name, self.readNoise))
    
    def initialCamera(self,cmd):
        """ Initial the MCS camera. """

        cmd.inform('text="Starting camera initialization."')
        p = sub.Popen(['/opt/EDTpdv/initcam', '-f', '/home/chyan/Canon50M/canon-8960x5778.cfg'],stdout=sub.PIPE,stderr=sub.PIPE)
        output, errors = p.communicate()
        string=errors[23:-1]
        if (string == 'done'):
            cmd.inform('text="Camera initialization message: %s"' % (string))
    
    def setExposureTime(self, cmd, expTime):
        """ Initial the MCS camera. """
        self.expTime = expTime
        cmd.inform('expTime=%f ms' % (expTime))
        
    def expose(self, cmd, expTime, expType, filename):
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
            pass
            # The expTime unit is ms.
            #time.sleep((expTime/1000.0) + self._exposureOverheadTime())

        # Command camera to do exposure sequence
        slicename=filename[0:14]
        p = sub.Popen(['canonexp', '-f', slicename, '-t', str(expTime), '-c'],stdout=sub.PIPE,stderr=sub.PIPE)
        output, errors = p.communicate()
        if (output == 'done'):
            cmd.inform('exposureState="done"')       

        if cmd:
            cmd.inform('exposureState="reading"')

        f = pyfits.open('/home/pfs/mhs/devel/ics_mcsActor/coadd.fits')      
        image = f[0].data
        #image = numpy.random.normal(self.biasLevel, 
        #                            scale=self.readNoise, 
        #                            size=self.imageSize).astype('u2')

        if expType != 'test':
            time.sleep(self._readoutTime())
        return image
        
        
