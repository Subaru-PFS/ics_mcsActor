#!/usr/bin/env python

from __future__ import print_function
from builtins import zip

from builtins import range
from builtins import object
import matplotlib
import matplotlib.pyplot as plt
import time

import os
import base64
import numpy
import astropy.io.fits as pyfits
import sys

import opscore.protocols.keys as keys
import opscore.protocols.types as types

from opscore.utility.qstr import qstr

import psycopg2
import psycopg2.extras
from xml.etree.ElementTree import dump

import mpfitCentroid.centroid as centroid

import numpy as np
import pylab as py

matplotlib.use('Agg')


class McsCmd(object):
    # Setting the default exposure time.
    dummy_filename = None

    def __init__(self, actor):
        # This lets us access the rest of the actor.
        self.actor = actor
        self.expTime = 1000

        self.simulationPath = None

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
            ('timeTest','',self.timeTest),
            ('seeingTest','',self.seeingTest),
            ('simulate', 'on <path>', self.simulateOn),
            ('simulate', 'off', self.simulateOff),
        ]

        # Define typed command arguments for the above commands.
        self.keys = keys.KeysDictionary("mcs_mcs", (1, 1),
                                        keys.Key("expTime", types.Float(), help="The exposure time"),
                                        keys.Key("expType", types.String(), help="The exposure type"),
                                        keys.Key("filename", types.String(), help="Image filename"),
                                        keys.Key("path", types.String(), help="Simulated image directory"),
                                        keys.Key("getArc", types.Int(), help="flag for arc image")
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

    def simulateOff(self, cmd):
        self.simulationPath = None
        cmd.finish('text="set simulation path to %s"' % str(self.simulationPath))

    def simulateOn(self, cmd):
        cmdKeys = cmd.cmd.keywords

        path = cmdKeys['path'].values[0]

        if not os.path.isdir(path):
            cmd.fail('text="path %s is not a directory"' % (path))
            return

        self.simulationPath = (path, 0)
        cmd.finish('text="set simulation path to %s"' % str(self.simulationPath))

    def getNextSimulationImage(self, cmd):
        import glob

        path, idx = self.simulationPath
        files = sorted(glob.glob(os.path.join(path, '*.fits')))
        if len(files) == 0:
            raise RuntimeError(f"no .fits files in {path}")

        if len(files) > idx+1:
            idx = 0

        imagePath = files[idx]
        image = pyfits.getdata(imagePath, 0)
        self.simulationPath = (path, idx+1)

        return image

    def getNextFilename(self, cmd):
        """ Fetch next image filename. 

        In real life, we will instantiate a Subaru-compliant image pathname generating object.  

        """

        ret = self.actor.cmdr.call(actor='gen2', cmdStr='getFrameId', timeLim=3.0)
        if ret.didFail:
            raise RuntimeError("getNextFilename failed!")

        self.actor.exposureID = self.actor.models['gen2'].keyVarDict['frameid'].valueList[0]

        path = os.path.join("$ICS_MHS_DATA_ROOT", 'mcs')
        path = os.path.expandvars(os.path.expanduser(path))

        if not os.path.isdir(path):
            os.makedirs(path, 0o755)
            
        return os.path.join(path, 'PFSC%06d00.fits' % (self.actor.exposureID))

    def _constructHeader(self, expType, expTime):
        ret = self.actor.cmdr.call(actor='gen2',
                                   cmdStr=f'getFitsCards \
                                            frameid={self.actor.exposureID} \
                                            expType={expType} expTime={expTime}',
                                   timeLim=3.0)
        if ret.didFail:
            raise RuntimeError("getFitsCards failed!")

        hdr = self.actor.models['gen2'].keyVarDict['header'].valueList[0]

        return pyfits.Header.fromstring(hdr)

    def getNextDummyFilename(self, cmd):
        """ Fetch next image filename. 

        In real life, we will instantiate a Subaru-compliant image pathname generating object.  

        """

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

        if self.simulationPath is None:
            image = self.actor.camera.expose(cmd, expTime, expType)
        else:
            image = self.getNextSimulationImage(cmd)

        hdr = self._constructHeader(expType, expTime)
        phdu = pyfits.PrimaryHDU(header=hdr)
        imgHdu = pyfits.CompImageHDU(image, compression_type='GZIP_2')
        hduList = pyfits.HDUList([phdu, imgHdu])

        hduList.writeto(filename, checksum=False, overwrite=True)

        cmd.inform('filename="%s"' % (filename))

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


    def doCentroidCoarse(self, cmd):
        

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

        
        cmd.finish('Statistics Calculated')
        
    def quickPlot(self,cmd):
        py.clf()
        npoint=len(self.actor.homes)//2
        for i in range(0,npoint):
            py.plot([self.actor.homes[i]],self.actor.homes[i+npoint],'dg')
        py.title("Centroids")

        py.savefig("test1.jpg")

    def fakeCentroidOnly(self,cmd):

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
        print(centroid.__file__)

        #cmd.inform('state="taking exposure"')
                
        #filename, image = self._doExpose(cmd, expTime, expType)
        
        #image=self._doFakeExpose(cmd, expTime, expType, "/Users/karr/GoogleDrive/TestData/home",0)
        
        # The encoding scheme is temporary, and will become encapsulated.
        cmd.inform('state="measuring"')

        #centroids = numpy.random.random(4800).astype('f4').reshape(2400,2)
        #self.dumpCentroidtoDB(cmd, centroids)

        cmd.inform('text="size = %s." '% (type(self.actor.image.astype('<i4'))))

        a = centroid.get_homes_call(self.actor.image.astype('<i4'))
        
        #a=get_homes_call(self.actor.image.astype('<i4'))
        homes=np.frombuffer(a,dtype='<f8')

        #centroidsStr = self._encodeArray(centroids)
        #cmd.inform('state="measured"; centroidsChunk=%s' % (centroidsStr))
        #
        cmd.inform('text="size = %d." '% (homes.shape))
        cmd.inform('state="centroids measured"')
        self.actor.homes=homes
        npoint=len(self.actor.homes)//2
        for i in range(0,npoint):
            print(self.actor.homes[i],self.actor.homes[i+npoint],'dg')
            cmd.inform('text="size = %f %f." '% (self.actor.homes[i],self.actor.homes[i+npoint]))

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
        
        image=self._doFakeExpose(cmd, expTime, expType, "/Users/karr/GoogleDrive/TestData/home",0)

        print("Read Image\n");

        #centroid call
        
        a = centroid.get_homes_call(image)

        #convert cython output into numpy array
        
        homes=np.frombuffer(a,dtype='<f8')

        print(homes)

        #second step: short exposure with arc image
        #image=pyfits.getdata('/Users/karr/GoogleDrive/first_move.fits').astype('<i4')
        #arc_image=pyfits.getdata('/Users/karr/GoogleDrive/first_move_arc.fits').astype('<i4')

        expTime=0.5
        image, arc_image=self._doFakeExpose(cmd, expTime, expType, "/Users/karr/GoogleDrive/first_move",1)
        
        #Call the centroiding/finding
        
        b = centroid.centroid_coarse_call(image,arc_image,homes)

        #convert from cython output to numpy typed array
        
        homepos=np.frombuffer(b,dtype=[('xp','<f8'),('yp','<f8'),('xt','<f8'),('yt','<f8'),('xc','<f8'),('yc','<f8'),('x','<f8'),('y','<f8'),('peak','<f8'),('back','<f8'),('fx','<f8'),('fy','<f8'),('qual','<f4'),('idnum','<f4')])

        #second move, same as the first

        #image=pyfits.getdata('./second_move.fits').astype('<i4')
        #arc_image=pyfits.getdata('./second_move.fits').astype('<i4')

        expTime=0.5
        image, arc_image=self._doFakeExpose(cmd, expTime, expType, "/Users/karr/GoogleDrive/second_move",1)

        b = centroid.centroid_coarse_call(image,arc_image,homes)

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
        c = centroid.centroid_fine_call(image,homes,xp,yp)

        #cython to numpy
        
        homepos=np.frombuffer(c,dtype=[('xp','<f8'),('yp','<f8'),('xt','<f8'),('yt','<f8'),('xc','<f8'),('yc','<f8'),('x','<f8'),('y','<f8'),('peak','<f8'),('back','<f8'),('fx','<f8'),('fy','<f8'),('qual','<f4'),('idnum','<f4')])

    def _makeTables(self, conn, doDrop=False):
        """ Create test tables. Someone said 20 measured things. """
        
        cmd = '''create table mcsPerFiber (
        id SERIAL PRIMARY KEY,
        frameId integer,
        moveId smallint,
        fiberId smallint,
        
        centroidX real, centroidY real,
        f1x real, f1y real,
        f2x real, f2y real,
        f3x real, f3y real,
        f4x real, f4y real,
        f5x real, f5y real,
        f6x real, f6y real,
        f7x real, f7y real,
        f8x real, f8y real,
        f9x real, f9y real
        );'''
    
        
        with conn.cursor() as curs:
            if doDrop:
                curs.execute('drop table mcsPerFiber')
            curs.execute(cmd)
        conn.commit()
    def _writeCentroids(self, centArr, nextRowId, frameId, moveId, conn=None):
        """ Write all measurements for a given (frameId, moveId) """
    
        # Save measurements to a CSV buffer
        measBuf = io.StringIO()
        np.savetxt(measBuf, centArr, delimiter=',', fmt='%0.6g')
        measBuf.seek(0,0)
        
        buf = io.StringIO()
        for l_i in range(len(centArr)):
            line = '%d,%d,%d,%d,%s' % (nextRowId + l_i, frameId, moveId, l_i,
                                       measBuf.readline())
            buf.write(line)
        buf.seek(0,0)
        
        if conn is not None:
            with conn.cursor() as curs:
                curs.copy_from(buf,'mcsPerFiber',',')
            conn.commit()
            buf.seek(0,0)
            
        return buf

    def _readCentroids(self, conn, frameId, moveId):
        """ Read all measurements for a given (frameId, moveId)"""
        buf = io.StringIO()
    
        cmd = f"""copy (select * from mcsPerFiber where frameId={frameId} and moveId={moveId}) to stdout delimiter ',' """
        with conn.cursor() as curs:
            curs.copy_expert(cmd, buf)
        conn.commit()
        buf.seek(0,0)
    
        # Skip the frameId, etc. columns.
        arr = np.genfromtxt(buf, dtype='f4',
                            delimiter=',',usecols=range(4,24))
        return arr

    def _writeDBtable(self):
        conn = psycopg2.connect("dbname='fps' user='pfs' host='localhost' password='pfspass'")
        cur = conn.cursor()
        cur.execute("select * from information_schema.tables where table_name=%s", ('test',))
        print(bool(cur.rowcount))
        conn.close()
        
    def timeTest(self,cmd):


        fwhm=3.
        hmin=3000
        boxsize=9
        expTime=1000.
        expType="object"
        filename, image = self._doExpose(cmd, expTime, expType)
        self.actor.image = image.astype('<i4')

        t1=time.time()
        for i in range(20):
            #expTime=1000.
            #expType="object"
            #filename, image = self._doExpose(cmd, expTime, expType)
            #self.actor.image = image.astype('<i4')

            a = centroid.centroid_only(self.actor.image,fwhm,hmin,boxsize)
        
            centroids=np.frombuffer(a,dtype='<f8')
            numpoints=len(centroids)//7
            self.actor.centroids=np.reshape(centroids,(numpoints,7))
            cmd.inform('text="size = %d." '% (numpoints))
        t2=time.time()
        cmd.inform('text="time = %f." '% ((t2-t1)/20.))

        cmd.finish('exposureState=done')
        
        
    def seeingTest(self,cmd):

        fwhm=3.
        hmin=3000
        boxsize=9
        expType="object"


        exptimes=np.array([500,1000.,2000.,3000.,4000.,5000,10000.])
        for expTime in exptimes:

            
            for i in range(30):
            
                filename, image = self._doExpose(cmd, expTime, expType)
                if(i==0):
                    cmd.inform('="expTime = %f. first= %s" '% (expTime,filename))
                    #self.actor.image = image.astype('<i4')
            cmd.inform('="expTime = %f. last= %s" '% (expTime,filename))

        cmd.finish('exposureState=done')
