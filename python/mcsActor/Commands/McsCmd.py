#!/usr/bin/env python

import os
import base64
import numpy
import astropy.io.fits as pyfits

import opscore.protocols.keys as keys
import opscore.protocols.types as types
#
from opscore.utility.qstr import qstr

#For Centroiding
from centroid import get_homes_call
from centroid import centroid_coarse_call
from centroid import centroid_fine_call
import pyfits
import numpy as np
import pylab as py


class McsCmd(object):
    # Setting the default exposure time.
    

    def __init__(self, actor):
        # This lets us access the rest of the actor.
        self.actor = actor
        self.expTime = 800

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
            ('expose_standard', '', self.expose),
            ('expose_long', '', self.expose),
            ('centroid', '<expTime>', self.centroid),
            ('test_centroid', '<expTime>', self.centroid),
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

        cmd.finish('text="camera connected!"')
        
    def status(self, cmd):
        """Report status and version; obtain and send current data"""

        self.actor.sendVersionKey(cmd)
        self.actor.camera.sendStatusKeys(cmd)
        self.actor.connectCamera(cmd)
        self.actor.camera.setExposureTime(cmd,self.expTime)

        cmd.inform('text="Present!"')
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

    def _doFakeExpose(self, cmd, expTime, expType, filename, getArc):
        
        """ Fake exposure, returns either image, or image + arc image. """

        #read the image file
        
        image=pyfits.getdata(filename+".fits",1)

        #if needed, read the arc file
        
        if(getArc==1):
            arc_image=pyfits.getdata(filename+"_arc.fits",1)
            return image, arc_image
        else:
            return image
        
    def _encodeArray(self, array):
        """ Temporarily wrap a binary array transfer encoding. """

        # Actually, we want dtype,naxis,axNlen,base64(array)
        return base64.b64encode(array.tostring())
        
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

    def test_centroid(self, cmd):
        
        """ 

        Demo Command to run a centroid sequence. 

        This needs to be split into a series of commands, with input from FPS,
        and appropriate handling of intput/output/configuration by either keywords
        or database, as decided. 

        """

        #Read in Simulated Data

        #First step: get the centroids in the home position, from a
        #long exposure. This does not do fibre identification

        #Fake camera image
        #image=pyfits.getdata('TestData/home.fits').astype('<i4')

        #get an image from a file, no arc image, long exposure time. 
        expTime=1
        expType='object'
        
        image=_doFakeExpose(cmd, expTime, expType, "TestData/home",0)

        #centroid call
        
        a=get_homes_call(image)

        #convert cython output into numpy array
        
        homes=np.frombuffer(a,dtype='<f8')

        #second step: short exposure with arc image
        #image=pyfits.getdata('TestData/first_move.fits').astype('<i4')
        #arc_image=pyfits.getdata('TestData/first_move_arc.fits').astype('<i4')

        expTime=0.5
        image, arc_image=_doFakeExpose(cmd, expTime, expType, "TestData/first_move",1)
        
        #Call the centroiding/finding
        
        b=centroid_coarse_call(image,arc_image,homes)

        #convert from cython output to numpy typed array
        
        homepos=np.frombuffer(b,dtype=[('xp','<f8'),('yp','<f8'),('xt','<f8'),('yt','<f8'),('xc','<f8'),('yc','<f8'),('x','<f8'),('y','<f8'),('peak','<f8'),('back','<f8'),('fx','<f8'),('fy','<f8'),('qual','<f4'),('idnum','<f4')])

        #second move, same as the first

        #image=pyfits.getdata('./second_move.fits').astype('<i4')
        #arc_image=pyfits.getdata('./second_move.fits').astype('<i4')

        expTime=0.5
        image, arc_image=_doFakeExpose(cmd, expTime, expType, "TestData/second_move",1)

        b=centroid_coarse_call(image,arc_image,homes)

        homepos=np.frombuffer(b,dtype=[('xp','<f8'),('yp','<f8'),('xt','<f8'),('yt','<f8'),('xc','<f8'),('yc','<f8'),('x','<f8'),('y','<f8'),('peak','<f8'),('back','<f8'),('fx','<f8'),('fy','<f8'),('qual','<f4'),('idnum','<f4')])


        #now for a fine move: long exposure, not arc image

        image=pyfits.getdata('./third_move.fits').astype('<i4')

        expTime=1.0
        image = _doFakeExpose(cmd, expTime, expType, "TestData/third_move",0)

        #we need to pass it the list of previous positions as well
        
        npos=homepos.shape[0]

        xp=np.zeros((npos))
        yp=np.zeros((npos))

        for i in range(npos):
            xp[i]=homepos[i][6]
            yp[i]=homepos[i][7]

        #and the call
        c=centroid_fine_call(image,homes,xp,yp)

        #cython to numpy
        
        homepos=np.frombuffer(c,dtype=[('xp','<f8'),('yp','<f8'),('xt','<f8'),('yt','<f8'),('xc','<f8'),('yc','<f8'),('x','<f8'),('y','<f8'),('peak','<f8'),('back','<f8'),('fx','<f8'),('fy','<f8'),('qual','<f4'),('idnum','<f4')])

        #this section demonstrates a plot of the final positions, patrol regions and home
        #positions, with a line connecting each home/fibre

        #xc=np.zeros((npos))
        #yc=np.zeros((npos))
        #xh=np.zeros((npos))
        #yh=np.zeros((npos))
        #xp=np.zeros((npos))
        #yp=np.zeros((npos))
        #
        #print homepos[1]
        #
        #for i in range(npos):
        #    xp[i]=homepos[i][0]
        #    yp[i]=homepos[i][1]
        #    xc[i]=homepos[i][4]
        #    yc[i]=homepos[i][5]
        #    xh[i]=homepos[i][6]
        #    yh[i]=homepos[i][7]
        #
        #
        #py.clf()
        #py.plot(xc,yc,'dg')
        #py.plot(xh,yh,'ob')
        #fig=py.gcf()
        #ax=fig.gca()
        #
        #
        #for i in range(len(xc)):
        #    circle=py.Circle((xh[i],yh[i]),55,color="black",fill=False)
        #    ax.add_artist(circle)
        #
        #    py.plot([xh[i],xc[i]],[yh[i],yc[i]])
        #
        #
        #py.show()
