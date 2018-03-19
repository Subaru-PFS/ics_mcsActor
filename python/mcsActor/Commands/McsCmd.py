#!/usr/bin/env python



from __future__ import print_function
from builtins import zip

from builtins import range
from builtins import object
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('Agg')

import os
import base64
import numpy
import astropy.io.fits as pyfits
import sys
import time

import opscore.protocols.keys as keys
import opscore.protocols.types as types

from opscore.utility.qstr import qstr

import psycopg2
import psycopg2.extras
from xml.etree.ElementTree import dump


sys.path.append("/home/pfs/mhs/devel/ics_mcsActor/python/mcsActor/mpfitCentroid")
from centroid import get_homes_call
from centroid import centroid_coarse_call
from centroid import centroid_fine_call

import pyfits
import numpy as np
import pylab as py
import centroid

class McsCmd(object):
    # Setting the default exposure time.
    dummy_filename = None

    def __init__(self, actor):
        # This lets us access the rest of the actor.
        self.actor = actor
        self.expTime = 1000

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
            ('centroidOnly', '<expTime>', self.centroidOnly),
            ('fakeCentroidOnly', '<expTime>', self.fakeCentroidOnly),
            ('test_centroid', '', self.test_centroid),
            ('fakeExpose','<expTime> <expType> <filename> <getArc>', self.fakeExpose),
            ('reconnect', '', self.reconnect),
            ('imageStats', '', self.imageStats),
            ('quickPlot', '', self.quickPlot),
        ]

        # Define typed command arguments for the above commands.
        self.keys = keys.KeysDictionary("mcs_mcs", (1, 1),
                                        keys.Key("expTime", types.Float(), help="The exposure time"),
                                        keys.Key("expType", types.String(), help="The exposure type"),
                                        keys.Key("filename", types.String(), help="Image filename"),
                                        keys.Key("getArc", types.Int(), help="flag for arc image")
                                        keys.Key("r1", types.Float(), help="lower bound for roundlim")
                                        keys.Key("r2", types.Float(), help="upper bound for roundlim")
                                        keys.Key("s1", types.Float(), help="lower bound for sharplim")
                                        keys.Key("s2", types.Float(), help="upper bound for sharplim")
                                        keys.Key("boxsize", types.Int(), help="fitting boxsize")
                                        keys.Key("fwhm", types.Float(), help="expected fwhm")
                                        keys.Key("thresh", types.Int(), help="threshold for peak finding")
                                        )

    def ping(self, cmd):
        """Query the actor for liveness/happiness."""
        self.actor.connectCamera(cmd)
        self.actor.camera.setExposureTime(cmd,self.expTime)

        cmd.finish("text='Present and (probably) well'")

    def reconnect(self, cmd):
        self.actor.connectCamera(cmd)
        self.actor.camera.setExposureTime(cmd,self.expTime)

        cmd.finish('text="Camera connected!"')
        
    def status(self, cmd):
        """Report status and version; obtain and send current data"""

        self.actor.sendVersionKey(cmd)
        self.actor.camera.sendStatusKeys(cmd)
        self.actor.connectCamera(cmd)
        self.actor.camera.setExposureTime(cmd,self.expTime)

        cmd.inform('text="MCS camera present!"')
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
    
    def dumpCentroidtoDB(self, cmd):
        """Connect to database and return json string to an attribute."""
        file = open("/home/pfs/mhs/devel/ics_mcsActor/etc/dbpasswd.cfg", "r")
        passstring = file.read() 
        cmd.inform('text="Connected to FPS database with pw %s."'%(passstring))
        try:
            conn = psycopg2.connect("dbname='fps' user='pfs' host='localhost' password="+passstring)
            cmd.diag('text="Connected to FPS database."')
        except:
            cmd.diag('text="I am unable to connect to the database."')
    #        print("I am unable to connect to the database.")
        #pass        
        cur = conn.cursor()
    
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
        self.actor.image = image
        
        #import pdb; pdb.set_trace()
        #plt.ion()
        #plt.hist(self.actor.image.ravel())
        #plt.savefig('foo.pdf')
        #plt.show()
        basename=filename[0:37]
        self.imageStats(cmd, basename)
        
        self.dumpCentroidtoDB(cmd)
        cmd.finish('exposureState=done')

    def manualInitCentroid(self,cmd):

        #a routine to manually adjust the values set in initCentroid,
        #for comisssioning purposes

        self.actor.fwhm = cmd.cmd.keywords["fwhm"].values[0]
        self.actor.fwhm = cmd.cmd.keywords["thresh"].values[0]
        self.actor.fwhm = cmd.cmd.keywords["boxsize"].values[0]
        r1 = cmd.cmd.keywords["r1"].values[0]
        r2 = cmd.cmd.keywords["f2"].values[0]
        s1 = cmd.cmd.keywords["s1"].values[0]
        s2 = cmd.cmd.keywords["s2"].values[0]

        self.actor.roundlim=[r1,r2]
        self.actor.sharplim=[s1,s2]
        
    def initCentroid(self,cmd):

        #initialize centroid parameters
        
        self.actor.fwhm=3.
        self.actor.thresh=1800
        self.actor.boxsize=9
        self.actor.roundlim=[-1.,1.]
        self.actor.sharplim=[0.05,0.5]

    def doCentroidCoarse(self, cmd):
        
        expTime = cmd.cmd.keywords["expTime"].values[0]
        expType = 'object' 

        
        
        pass
        
    def _doCentroid(self,cmd,image,fittype):

        pass

        
    def fakeExpose(self,cmd):
        
        cmdKeys = cmd.cmd.keywords

        expTime = cmdKeys['expTime'].values[0]
        expType = cmdKeys['expType'].values[0]
        filename = cmdKeys['filename'].values[0]
        getArc = cmdKeys['getArc'].values[0]

        if(getArc==1):
            image, arc_image=self._doFakeExpose(cmd,expTime,expType,filename,getArc)
        else:
            image = self._doFakeExpose(cmd,expTime,expType,filename,getArc)

        self.actor.image = image

        cmd.finish('exposureState=done')
            
    def _doFakeExpose(self, cmd,expTime,expType,filename,getArc):
        
        
        """ Fake exposure, returns either image, or image + arc image. """

        # Actually, we want dtype,naxis,axNlen,base64(array)
        return base64.b64encode(array.tostring())
        

        image=pyfits.getdata(filename+".fits",0).astype('<i4')
        
        #if needed, read the arc file
        
        if(getArc==1):
            arc_image=pyfits.getdata(filename+"_arc.fits",0).astype('<i4')
            return image, arc_image
        else:
            return image
 
    def _encodeArray(self, array):
        """ Temporarily wrap a binary array transfer encoding. """

        # Actually, we want dtype,naxis,axNlen,base64(array)
        return base64.b64encode(array.tostring())

    def imageStats(self, cmd, basename):

        cmd.inform('text="image median = %d." '% (np.median(self.actor.image))) 
        cmd.inform('text="image mean = %d." '% (self.actor.image.mean())) 
        cmd.inform('text="image min = %d." '% (self.actor.image.min())) 
        cmd.inform('text="image max = %d." '% (self.actor.image.max())) 
        
        #plt.hist(self.actor.image.ravel())
        #plt.savefig(basename+'.pdf')
        #plt.show()

        
        cmd.finish('Statistics Calculated')
        
    def quickPlot(self,cmd):
        py.clf()
        npoint=len(self.actor.homes)//2
        for i in range(0,npoint):
            py.plot([self.actor.homes[i]],self.actor.homes[i+npoint],'dg')
        py.title("Centroids")

        py.savefig("test1.jpg")

def fakeCentroidOnly(self,cmd):

        #generates test data of the right format, useful for
        #testing database injection
    
        cmd.inform('state="measuring"')

        npos=2350
        centroids=np.zeros((2350,7))

        #fake positions
        pos=np.meshgrid(np.arange(50),np.arange(47))
        centroids[:,0]=pos[0]*150
        centroids[:,1]=pos[1]*100

        #fake FWHM
        centroids[:,2]=np.normal(3.,0.4,npos)
        centroids[:,3]=np.normal(3.,0.4,npos)

        #fake peaks
        centroids[:,4]=np.normal(5000,100,npos)

        #fake backgrounds
        centroids[:,5]=np.normal(800,30,npos)

        #fake qualities
        centroids[:,6]=np.ones(npos)

        self.actor.centroids=centroids

        cmd.inform('state="finished"')

    def centroidOnly(self, cmd):
        """ Take an exposure and measure centroids. """
        
        expTime = cmd.cmd.keywords["expTime"].values[0]
        expType = 'object' 
        
        cmd.inform('state="measuring"')

        #call teh centroiding routine. The parameters need to have been set with initCentroid
        #or manualInitCentroid prior to this. 
        
        a=centroid_only(self.actor.image.astype('<i4'),self.actor.fwhm,self.actor.thresh,self.actor.boxsize,self.actor.roundlim,self.actor.sharplim)

        #read into a numpy buffer
        centroids=np.frombuffer(a,dtype='<f8')

        #reshape
        np=len(centroids)//7
        self.actor.centroids=np.reshape(centroids,(np,7))

        #a=get_homes_call(self.actor.image.astype('<i4'))
        #homes=np.frombuffer(a,dtype='<f8')

        #update status
        cmd.inform('state="centroids measured"')
        cmd.inform('text="detected %d points" '% (np))
        
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
        
        image=self._doFakeExpose(cmd, expTime, expType, "/Users/karr/GoogleDrive/TestData/home",0)

        print("Read Image\n");

        #centroid call
        
        a=get_homes_call(image)

        #convert cython output into numpy array
        
        homes=np.frombuffer(a,dtype='<f8')

        print(homes)

        #second step: short exposure with arc image
        #image=pyfits.getdata('/Users/karr/GoogleDrive/first_move.fits').astype('<i4')
        #arc_image=pyfits.getdata('/Users/karr/GoogleDrive/first_move_arc.fits').astype('<i4')

        expTime=0.5
        image, arc_image=self._doFakeExpose(cmd, expTime, expType, "/Users/karr/GoogleDrive/first_move",1)
        
        filename, image, arcimage = self._doExpose(cmd, expTime, expType)

        #Call the centroiding/finding
        
        b=centroid_coarse_call(image,arc_image,homes)

        #convert from cython output to numpy typed array
        
        fibreid=np.frombuffer(b,dtype=[('xp','<f8'),('yp','<f8'),('xt','<f8'),('yt','<f8'),('xc','<f8'),('yc','<f8'),('x','<f8'),('y','<f8'),('peak','<f8'),('back','<f8'),('fx','<f8'),('fy','<f8'),('qual','<f4'),('idnum','<f4')])

        #second move, same as the first

        #image=pyfits.getdata('./second_move.fits').astype('<i4')
        #arc_image=pyfits.getdata('./second_move.fits').astype('<i4')

        expTime=0.5
        image, arc_image=self._doFakeExpose(cmd, expTime, expType, "/Users/karr/GoogleDrive/second_move",1)

        b=centroid_coarse_call(image,arc_image,homes)

        homepos=np.frombuffer(b,dtype=[('xp','<f8'),('yp','<f8'),('xt','<f8'),('yt','<f8'),('xc','<f8'),('yc','<f8'),('x','<f8'),('y','<f8'),('peak','<f8'),('back','<f8'),('fx','<f8'),('fy','<f8'),('qual','<f4'),('idnum','<f4')])


        #now for a fine move: long exposure, not arc image

        image=pyfits.getdata('./third_move.fits').astype('<i4')

        expTime=1.0
        image = _doFakeExpose(cmd, expTime, expType, "/Users/karr/GoogleDrive/third_move",0)

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

    def displayImageResults(self, cmd):

        """

        A simple routine to plot the image with the detected centroids
        overlaid, for calibration purposes.

        """
        
        filename = cmd.cmd.keywords["fileName"].values[0]
        image = pyfits.getdata(filename)

        py.imshow(image)
        py.scatter(self.actor.centroids[:,0],self.actor.centroids[:,1],'og')
        py.show()

    def showCentroidStats():

        """

        a simple routine to display various parameters from a centroid
        run, primarily for calibration/testing purposes

        """

        av_rms_x=self.actor.centroids[:,2].mean()
        av_rms_y=self.actor.centroids[:,3].mean()
        av_peak=self.actor.centroids[:,4].mean()
        av_back=self.actor.centroids[:,5].mean()

        
    def rmsValues():

        """

        wrapper routine to pull the centroids of an image sequence
        from the database and calculate RMS values of the position,
        along with some other diagnostic values.

        """

        centroidList=[]
        n=10
        
        for i in range(n):
            
            filename, image = self._doExpose(cmd, expTime, expType)
            self.actor.image = image
        
            a=centroid_only(self.actor.image.astype('<i4'),self.actor.fwhm,self.actor.thresh,self.actor.boxsize,self.actor.roundlim,self.actor.sharplim)
            np=len(centroids)//7
            centroids=np.reshape(centroids,(np,7))

            centroidList.append(centroids)

        rmsPlots(centroidList)

    _doFlatField(self,cmd,image,flat):
        
        
    def timeTest(self, cmd):

        time1=time.time()
        for i in range(20):
            filename, image = self._doExpose(cmd, expTime, expType)
            self.actor.image = image
        
            a=centroid_only(self.actor.image.astype('<i4'),self.actor.fwhm,self.actor.thresh,self.actor.boxsize,self.actor.roundlim,self.actor.sharplim)
        
            #read into a numpy buffer
            centroids=np.frombuffer(a,dtype='<f8')

            #reshape
            np=len(centroids)//7
            self.actor.centroids=np.reshape(centroids,(np,7))
        time2=time.time()
