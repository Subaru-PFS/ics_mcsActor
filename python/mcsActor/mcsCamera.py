from builtins import object
import numpy
import time
import subprocess as sub
import astropy.io.fits as pyfits
import os
import shutil


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

        cmd.inform('text="Starting camera initialization."')
        p = sub.Popen(['/opt/EDTpdv/initcam', '-f',  configPath],
                      stdout=sub.PIPE, stderr=sub.PIPE)
        output, errors = p.communicate()
        string = errors[23:-1]
        if (string == 'done'):
            cmd.inform('text="Camera initialization message: %s"' % (string))

    def setExposureTime(self, cmd, expTime):
        """ Initial the MCS camera. """
        self.expTime = expTime
        cmd.inform('expTime=%f ms' % (expTime))

    def expose(self, cmd, expTime, expType, filename, doCopy=True, flip=False):
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

        if not expType:
            expType = 'test'
        if cmd:
            cmd.inform('exposureState="exposing"')

        if expType in ('dark', 'bias') and expTime > 0:
            cmd.inform('text="ETYPE: %s"' % (expType))
            t1 = time.time()

            # Command camera to do exposure sequence
            slicename = filename[0:33]+'_'
            cmd.inform('text="slice name: %s"' % (slicename))
            p = sub.Popen(['canonexp', '-e', expType, '-f', slicename, '-t', str(expTime), '-c'],
                          cwd=self.coaddDir, bufsize=1, stdout=sub.PIPE, stderr=sub.PIPE)
            output, errors = p.communicate()
            t2 = time.time()

        if expType in ('flat') and expTime > 0:
            cmd.inform('text="ETYPE: %s"' % (expType))
            t1 = time.time()
            slicename = filename[0:33]+'_'
            cmd.inform('text="slice name: %s"' % (slicename))
            cmd.inform('text="open shutter"')
            p = sub.Popen(['shutter', '-o'], bufsize=1, stdout=sub.PIPE, stderr=sub.PIPE)
            output, errors = p.communicate()

            expType = 'flat'
            p = sub.Popen(['canonexp', '-e', expType, '-f', slicename, '-t', str(expTime), '-c'],
                          cmd=self.coaddDir, bufsize=1, stdout=sub.PIPE, stderr=sub.PIPE)
            output, errors = p.communicate()
            cmd.inform('text="taking exposure"')

            p = sub.Popen(['shutter', '-c'], bufsize=1, stdout=sub.PIPE, stderr=sub.PIPE)
            output, errors = p.communicate()
            cmd.inform('text="shutter closed"')
            t2 = time.time()

        if expType in ('test', 'object') and expTime > 0:
            t1 = time.time()

            # Command camera to do exposure sequence
            slicename = filename[0:33]+'_'
            cmd.inform('text="slice name: %s"' % (slicename))
            p = sub.Popen(['canonexp', '-f', slicename, '-t', str(expTime), '-c'],
                          cwd=self.coaddDir, bufsize=1, stdout=sub.PIPE, stderr=sub.PIPE)
            output, errors = p.communicate()
            t2 = time.time()

        if (output == 'Exposure sequence is done'):
            cmd.inform('exposureState="done"')

        if cmd:
            cmd.inform('exposureState="reading"')

        coaddpath = os.path.join(self.coaddDir, 'coadd.fits')
        if doCopy:
            try:
                shutil.copy(coaddpath, filename)
            except Exception as e:
                cmd.warn(f'text="failed to copy coadd {coaddpath} to {filename}: {e}"')

        f = pyfits.open(coaddpath)

        # Remove the temporary file
        try:
            os.remove(coaddpath)
        except OSError:
            pass

        image = f[0].data
        t3 = time.time()
        cmd.inform('text="Time for exposure = %f." ' % ((t2-t1)/1.))
        cmd.inform('text="Time for image loading= %f." ' % ((t3-t2)/1.))

        if flip is True:
            image=(np.flip(image).T).copy(order='C')

        return image
