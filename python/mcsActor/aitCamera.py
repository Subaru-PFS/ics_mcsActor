from __future__ import division
from builtins import str
from builtins import object
import subprocess as sub
import numpy
import time
import astropy.io.fits as pyfits


class Camera(object):
    pass


class aitCamera(Camera):
    def __init__(self):
        # Should read from config file....
        self.imageSize = (10000, 7096)
        self.biasLevel = 100
        self.readNoise = 5.0
        self.name = 'aitCamera'

    def _readoutTime(self):
        return 0.5

    def _exposureOverheadTime(self):
        return 0.1

    def sendStatusKeys(self, cmd):
        """ Send our status keys to the given command. """

        cmd.inform('cameraName=%s; readNoise=%0.2f' % (self.name, self.readNoise))

    def initialCamera(self, cmd):
        """ Initial the AIT camera. """

        cmd.inform('text="Starting camera initialization."')
        p = sub.Popen(['/opt/EDTpdv/initcam', '-f', '/home/chyan/71MP/illunis-71mp.cfg'],
                      stdout=sub.PIPE, stderr=sub.PIPE)
        output, errors = p.communicate()
        string = errors[23:-1]
        if (string == 'done'):
            cmd.inform('text="Camera initialization message: %s"' % (string))

    def setExposureTime(self, cmd, expTime):
        """ Initial the AIT camera. """
        p = sub.Popen(['rmodcontrol', '-e', str(expTime)], stdout=sub.PIPE, stderr=sub.PIPE)
        output, errors = p.communicate()
        if (str(expTime) == output[28:-6]):
            cmd.inform('expTime=%f ms' % (int(output[28:-6])))

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
            # The expTime unit is ms.
            time.sleep((expTime / 1000.0) + self._exposureOverheadTime())

        # Command camera to do exposure sequence.
        p = sub.Popen(['rmodexposure', '-f', filename], stdout=sub.PIPE, stderr=sub.PIPE)
        output, errors = p.communicate()
        if (output == 'done'):
            cmd.inform('exposureState="done"')

        if cmd:
            cmd.inform('exposureState="reading"')

        f = pyfits.open('/home/chyan/mhs/data/mcs/schmidt_fiber_snr400_rmod71.fits')
        image = f[0].data
        # image = numpy.random.normal(self.biasLevel,
        #                            scale=self.readNoise,
        #                            size=self.imageSize).astype('u2')

        if expType != 'test':
            time.sleep(self._readoutTime())
        return image
