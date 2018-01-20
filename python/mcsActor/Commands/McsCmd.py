#!/usr/bin/env python

import os
import base64
import numpy
import astropy.io.fits as pyfits
import sys

import opscore.protocols.keys as keys
import opscore.protocols.types as types

from opscore.utility.qstr import qstr

sys.path.append("/home/chyan/mhs/devel/ics_mcsActor/python/mcsActor/mpfitCentroid")
from centroid import get_homes_call
from centroid import centroid_coarse_call
from centroid import centroid_fine_call


class McsCmd(object):
    # Setting the default exposure time.
    dummy_filename = None

    def __init__(self, actor):
        # This lets us access the rest of the actor.
        self.actor = actor
        self.expTime = 100

        # Declare the commands we implement. When the actor is started
        # these are registered with the parser, which will call the
        # associated methods when matched. The callbacks will be
        # passed a single argument, the parsed and typed command.
        #
        self.vocab = [
            ('ping', '', self.ping),
            ('status', '', self.status),
            ('mockexpose', '@(bias|test)', self.mockexpose),
            ('mockexpose', '@(dark|object) <expTime>', self.mockexpose),
            ('expose', '@(bias|test)', self.expose),
            ('expose', '@(dark|object) <expTime>', self.expose),
            ('centroid', '<expTime>', self.centroid),
            ('centroidOnDummy', '<expTime>', self.centroidOnDummy),
            ('reconnect', '', self.reconnect),
        ]

        # Define typed command arguments for the above commands.
        self.keys = keys.KeysDictionary("mcs_mcs", (1, 1),
                                        keys.Key("expTime", types.Float(), help="The exposure time"),
                                        )


    def ping(self, cmd):
        """Query the actor for liveness/happiness."""
        self.actor.connectCamera(cmd)
        self.actor.camera.setExposureTime(cmd,self.expTime)

        cmd.finish("text='Present and (probably) well'")

    def reconnect(self, cmd):
        self.actor.connectCamera(cmd)
        self.actor.camera.setExposureTime(cmd,self.expTime)

        cmd.finish('text="AIT camera connected!"')
        
    def status(self, cmd):
        """Report status and version; obtain and send current data"""

        self.actor.sendVersionKey(cmd)
        self.actor.camera.sendStatusKeys(cmd)
        self.actor.connectCamera(cmd)
        self.actor.camera.setExposureTime(cmd,self.expTime)

        cmd.inform('text="AIT camera present!"')
        cmd.finish()

    def getNextFilename(self, cmd):
        """ Fetch next image filename. 

        In real life, we will instantiate a Subaru-compliant image pathname generating object.  

        """
        
        self.actor.exposureID += 1
        path = os.path.join("$ICS_MHS_DATA_ROOT", 'mcs')
        path = os.path.expandvars(os.path.expanduser(path))

        if not os.path.isdir(path):
            os.makedirs(path, 0o755)
            
        return os.path.join(path, 'MCSA%010d.fits' % (self.actor.exposureID))

    def getNextDummyFilename(self, cmd):
        """ Fetch next image filename. 

        In real life, we will instantiate a Subaru-compliant image pathname generating object.  

        """
        
        #self.actor.exposureID += 1
        path = os.path.join("$ICS_MHS_DATA_ROOT", 'mcs')
        path = os.path.expandvars(os.path.expanduser(path))

        if not os.path.isdir(path):
            os.makedirs(path, 0o755)
            
        return os.path.join(path, 'dummy_MCSA%010d.fits' % (self.actor.exposureID))

    def _doMockExpose(self, cmd, expTime, expType):
        """ Take an exposure and save it to disk. """

        filename = self.getNextFilename(cmd)
        #dummy_filename = self.getNextDummyFilename(cmd)

        f = pyfits.open('/home/chyan/mhs/data/mcs/schmidt_fiber_snr400_rmod71.fits')      
        image = f[0].data
        #image = self.actor.camera.expose(cmd, expTime, expType, filename)
        pyfits.writeto(filename, image, checksum=False, clobber=True)
        cmd.inform("filename=%s " % (filename))

        return filename, image
   
    def _doExpose(self, cmd, expTime, expType):
        """ Take an exposure and save it to disk. """

        filename = self.getNextFilename(cmd)
        dummy_filename = self.getNextDummyFilename(cmd)
        self.dummy_filename = dummy_filename
        
        image = self.actor.camera.expose(cmd, expTime, expType, filename)
        pyfits.writeto(dummy_filename, image, checksum=False, clobber=True)
        cmd.inform("filename=%s and dummy file=%s" % (filename, dummy_filename))

        return filename, image
     
    def mockexpose(self, cmd):
        """ Take an exposure and return mock image. Does not centroid. """

        expType = cmd.cmd.keywords[0].name
        if expType in ('bias', 'test'):
            expTime = self.expTime
        else:
            expTime = cmd.cmd.keywords["expTime"].values[0]

        #if (expTime != self.expTime):
            #self.actor.camera.setExposureTime(cmd,expTime)
    

        cmd.diag('text="Exposure time now is %d ms." '% (expTime))    
        

        filename, image = self._doMockExpose(cmd, expTime, expType)
        cmd.finish('exposureState=done')
           
    def expose(self, cmd):
        """ Take an exposure. Does not centroid. """

        expType = cmd.cmd.keywords[0].name
        if expType in ('bias', 'test'):
            expTime = self.expTime
        else:
            expTime = cmd.cmd.keywords["expTime"].values[0]

        if (expTime != self.expTime):
            self.actor.camera.setExposureTime(cmd,expTime)
    

        cmd.diag('text="Exposure time now is %d ms." '% (expTime))    
        

        filename, image = self._doExpose(cmd, expTime, expType)
        cmd.finish('exposureState=done')

    def _encodeArray(self, array):
        """ Temporarily wrap a binary array transfer encoding. """

        # Actually, we want dtype,naxis,axNlen,base64(array)
        return base64.b64encode(array.tostring())
        
    def centroidOnDummy(self, cmd):
        """ Take an exposure and measure centroids. """

        expTime = cmd.cmd.keywords["expTime"].values[0]
        expType = 'object' if expTime > 0 else 'test'

        filename, image = self._doExpose(cmd, expTime, expType)

        # The encoding scheme is temporary, and will become encapsulated.
        cmd.inform('state="measuring"')
        
        image=pyfits.getdata(self.dummy_filename).astype('<i4')

        #get homes positions
        a=get_homes_call(image)

        #turn into numpy array
        home_centroids=np.frombuffer(a,dtype='<f8')
        
        #centroids = numpy.random.random(4800).astype('f4').reshape(2400,2)

        centroidsStr = self._encodeArray(home_centroids)
        cmd.inform('state="measured"; centroidsChunk=%s' % (centroidsStr))

        cmd.finish('exposureState=done')


    def centroid(self, cmd):
        """ Take an exposure and measure centroids. """

        expTime = cmd.cmd.keywords["expTime"].values[0]
        expType = 'object' if expTime > 0 else 'test'

        filename, image = self._doExpose(cmd, expTime, expType)

        # The encoding scheme is temporary, and will become encapsulated.
        cmd.inform('state="measuring"')
        centroids = numpy.random.random(4800).astype('f4').reshape(2400,2)

        centroidsStr = self._encodeArray(centroids)
        cmd.inform('state="measured"; centroidsChunk=%s' % (centroidsStr))

        cmd.finish('exposureState=done')

