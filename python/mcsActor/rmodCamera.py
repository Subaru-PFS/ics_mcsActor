from builtins import object
import numpy as np
import time
import subprocess as sub
import astropy.io.fits as pyfits
import os
import shutil


class Camera(object):
    pass


class rmodCamera(Camera):
    def __init__(self):
        # Should read from config file....
        self.imageSize = (8960, 5778)
        self.biasLevel = 100
        self.readNoise = 5.0
        self.name = 'RMOD_71M'
        self.expTimeDefault = 800  # ms
        self.expTime = 800  # ms
        self.coaddDir = '/tmp'

    def _readoutTime(self):
        return 0.5

    def _exposureOverheadTime(self):
        return 0.1

    def sendStatusKeys(self, cmd):
        """ Send our status keys to the given command. """

        cmd.inform('cameraName=%s; readNoise=%0.2f' % (self.name, self.readNoise))

    def initialCamera(self, cmd, configPath):
        """ Initial the MCS camera. """

        cmd.inform(f'text="Starting camera initialization from {configPath}"')
        p = sub.Popen(['/opt/EDTpdv/initcam', '-f', configPath],
                      stdout=sub.PIPE, stderr=sub.PIPE)
        output, errors = p.communicate()
        string = errors[23:-1]
        if (string == 'done'):
            cmd.inform('text="Camera initialization message: %s"' % (string))

    def setExposureTime(self, cmd, expTime):
        """ Initial the MCS camera. """
        self.expTime = expTime
        p = sub.Popen(['rmodcontrol', '-e', f'{expTime}'], stdout=sub.PIPE, stderr=sub.PIPE)
        output, errors = p.communicate()

        cmd.inform('expTime=%f ms' % (expTime))

    def expose(self, cmd, expTime, expType, filename, doCopy=True):
        """ Generate an 'exposure' image. We don't have an actual camera, so generate some 
        plausible image. 

        Args:
           cmd     - a Command object to report to. Ignored if None.
           expTime - the fake exposure time. 
           expType - ("bias", "dark", "object", "flat", "test")
           filename - the template filename
           doCopy  - if True, copy coadd into filename

        Returns:
           - the image.

        Keys:
           exposureState
        """
        self.setExposureTime(cmd, expTime)

        if not expType:
            expType = 'test'
        if cmd:
            cmd.inform('exposureState="exposing"')

        if expType in ('dark', 'bias') and expTime > 0:
            cmd.inform('text="ETYPE: %s"' % (expType))
            t1 = time.time()

            # Command camera to do exposure sequence
            slicename = filename[0:33]+'.fits'
            cmd.inform('text="slice name: %s"' % (slicename))
            p = sub.Popen(['rmodexposure', '-f', slicename, '-l', '1'], stdout=sub.PIPE, stderr=sub.PIPE)
            output, errors = p.communicate()
            t2 = time.time()

        if expType in ('test', 'object') and expTime > 0:
            t1 = time.time()

            # Command camera to do exposure sequence
            slicename = filename[0:33]+'_'
            cmd.inform('text="slice name: %s"' % (slicename))
            p = sub.Popen(['rmodexposure', '-f', slicename, '-l', '1'], stdout=sub.PIPE, stderr=sub.PIPE)

            output, errors = p.communicate()
            if errors is not None:
                cmd.warn(f'text="exposure command failed with error: {errors}"')

            t2 = time.time()

        if cmd:
            cmd.inform('exposureState="reading"')

        f = pyfits.open(slicename)

        image = f[0].data
        t3 = time.time()

        
        cmd.inform('text="Time for exposure = %f." ' % ((t2-t1)/1.))
        cmd.inform('text="Time for image loading= %f." ' % ((t3-t2)/1.))

        # Reset exposure time to default after readout.
        if self.expTime > self.expTimeDefault:
            self.setExposureTime(cmd, self.expTimeDefault)
            cmd.inform('text="Set to a default exposure time = %f." ' % (self.expTimeDefault))

        return image
