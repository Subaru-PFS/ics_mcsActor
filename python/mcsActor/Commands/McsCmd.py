#!/usr/bin/env python

from __future__ import print_function
from builtins import zip

from builtins import range
from builtins import object
import ast

import matplotlib
import matplotlib.pyplot as plt
import time
import datetime
import io
import queue
import threading

import os
import base64
import astropy.io.fits as pyfits
import fitsio
import sys

import opscore.protocols.keys as keys
import opscore.protocols.types as types

from actorcore.utility import fits as fitsUtils
from opscore.utility.qstr import qstr

import pfs.utils.coordinates.MakeWCSCards as pfsWcs
# import pfs.utils.config as pfsConfig

import psycopg2
import psycopg2.extras
from xml.etree.ElementTree import dump

try:
    import mcsActor.mpfitCentroid.centroid as centroid
except:
    pass

import numpy as np
import pylab as py

# matplotlib.use('Agg')
Bool=bool

class McsCmd(object):

    def __init__(self, actor):
        # This lets us access the rest of the actor.
        self.actor = actor
        self.expTime = 1000
        self.newTable = None
        self.simulationPath = None
        self.exposureID = None
        self._conn = None

        self.db='db-ics'
        self.setCentroidParams(None, doFinish=False)
        
        # Declare the commands we implement. When the actor is started
        # these are registered with the parser, which will call the
        # associated methods when matched. The callbacks will be
        # passed a single argument, the parsed and typed command.
        #
        self.vocab = [
            ('ping', '', self.ping),
            ('status', '', self.status),
            ('expose', '@(bias|test) [<visit>]', self.expose),
            ('expose', '@(dark|flat) <expTime> [<visit>]', self.expose),
            ('expose', '@object <expTime> [<visit>] [@doCentroid]', self.expose),
            ('runCentroid', '[@newTable]', self.runCentroid),
            ('createTables', '[@drop]', self.createTables),
            ('fakeCentroidOnly', '<expTime>', self.fakeCentroidOnly),
            ('test_centroid', '', self.test_centroid),
            ('reconnect', '', self.reconnect),
            ('imageStats', '', self.imageStats),
            ('quickPlot', '', self.quickPlot),
            ('timeTestFull','',self.timeTestFull),
            ('seeingTest','',self.seeingTest),
            ('setCentroidParams','[<fwhm>] [<thresh>] [<rl>] [<rh>] [<sl>] [<sh>] [<boxsize>]',
             self.setCentroidParams),
            ('simulate', '<path>', self.simulateOn),
            ('simulate', 'off', self.simulateOff),
        ]

        # Define typed command arguments for the above commands.
        self.keys = keys.KeysDictionary("mcs_mcs", (1, 1),
                                        keys.Key("expTime", types.Float(), help="The exposure time, seconds"),
                                        keys.Key("expType", types.String(), help="The exposure type"),
                                        keys.Key("visit", types.Int(), help="exposure visit"),
                                        keys.Key("path", types.String(), help="Simulated image directory"),
                                        keys.Key("getArc", types.Int(), help="flag for arc image"),
                                        keys.Key("fwhm", types.Float(), help="fwhm for centroid routine"),
                                        keys.Key("boxsize", types.Int(), help="boxsize for centroid routine"),
                                        keys.Key("thresh", types.Int(), help="thresh for centroid routine"),
                                        keys.Key("rl", types.Float(), help="rl for centroid routine"),
                                        keys.Key("rh", types.Float(), help="rh for centroid routine"),
                                        keys.Key("sl", types.Float(), help="sh for centroid routine"),
                                        keys.Key("sh", types.Float(), help="sh for centroid routine")
                                        )

    @property
    def conn(self):
        if self._conn is not None:
            return self._conn

        pwpath=os.path.join(os.environ['ICS_MCSACTOR_DIR'],
                            "etc", "dbpasswd.cfg")

        try:
            file = open(pwpath, "r")
            passstring = file.read()
        except:
            raise RuntimeError(f"could not get db password from {pwpath}")

        try:
            connString = "dbname='fps' user='pfs' host="+self.db+" password="+passstring
            self.actor.logger.info(f'connecting to {connString}')
            conn = psycopg2.connect(connString)
            self._conn = conn
        except Exception as e:
            raise RuntimeError("unable to connect to the database {connString}: {e}")

        return self._conn

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

        self.simulationPath = (path, 0, '')
        cmd.finish('text="set simulation path to %s"' % str(self.simulationPath))

    def getNextSimulationImage(self, cmd):
        import glob

        path, idx, lastname = self.simulationPath
        files = sorted(glob.glob(os.path.join(path, '*.fits')))
        cmd.debug('text="%i of %i files in %s..."' % (idx, len(files), path))
        if len(files) == 0:
            raise RuntimeError(f"no .fits files in {path}")

        if  idx+1 > len(files):
            idx = 0

        imagePath = files[idx]
        image = pyfits.getdata(imagePath, 0)
        self.simulationPath = (path, idx+1, imagePath)
        cmd.debug('text="returning simulation file %s"' % (imagePath))
        return image

    def requestNextFilename(self, cmd):
        """ Return a queue which will eventually contain a filename. """

        q = queue.Queue()

        def worker(q=q, cmd=cmd):
            filename = self.getNextFilename(cmd)
            q.put(filename)
            q.task_done()

        task = threading.Thread(target=worker,
                                name='visitFetcher', daemon=True)
        task.start()

        return q

    def getNextFilename(self, cmd):
        """ Fetch next image filename. 

        In real life, we will instantiate a Subaru-compliant image pathname generating object.  

        """

        if self.exposureID is None:
            ret = self.actor.cmdr.call(actor='gen2', cmdStr='getVisit', timeLim=10.0)
            if ret.didFail:
                raise RuntimeError("getNextFilename failed getting a visit number in 10s!")

            self.actor.exposureID = self.actor.models['gen2'].keyVarDict['visit'].valueList[0]

        path = os.path.join("$ICS_MHS_DATA_ROOT", 'mcs', time.strftime('%Y-%m-%d'))
        path = os.path.expandvars(os.path.expanduser(path))
        if not os.path.isdir(path):
            os.makedirs(path, 0o755)
            
        newpath = os.path.join(path, 'PFSC%06d00.fits' % (self.actor.exposureID))

        return newpath
    
    def _getInstHeader(self, cmd):
        """ Gather FITS cards from all actors we are interested in. """

        # For now, do _not_ add gen2 cards, since we still have the gen2Actor generate them.
        modelNames = set(self.actor.models.keys())
        modelNames.discard("gen2")
        
        cmd.debug('text="fetching MHS cards for %s..."' % (modelNames))
        cards = fitsUtils.gatherHeaderCards(cmd, self.actor,
                                            modelNames=modelNames, shortNames=True)
        cmd.debug('text="fetched %d MHS cards..."' % (len(cards)))

        # Until we convert to fitsio, convert cards to pyfits
        pycards = []
        for c in cards:
            if isinstance(c, str):
                pcard = 'COMMENT', c
            else:
                pcard = c['name'], c['value'], c.get('comment', '')
            pycards.append(pcard)
            cmd.debug('text=%s' % (qstr("fetched card: %s" % (str(pcard)))))

        if self.simulationPath is not None:
            pcard = 'W_MCSMNM', self.simulationPath[-1], 'Simulated file path'
            pycards.append(pcard)

        return pycards

    def _constructHeader(self, cmd, filename, expType, expTime):
        if expType == 'bias':
            expTime = 0.0

        # Subaru rules
        if expType == 'object':
            expType = 'acquisition'

        frameId = os.path.splitext(os.path.basename(filename))[0]
        ret = self.actor.cmdr.call(actor='gen2',
                                   cmdStr=f'getFitsCards \
                                            frameid={frameId} \
                                            expType={expType} expTime={expTime/1000.0}',
                                   timeLim=3.0)
        if ret.didFail:
            raise RuntimeError("getFitsCards failed!")

        hdrString = self.actor.models['gen2'].keyVarDict['header'].valueList[0]
        hdrString = base64.b64decode(hdrString).decode('latin-1')
        try:
            hdr = pyfits.Header.fromstring(hdrString)
        except Exception as e:
            cmd.warn('text="FAILED to fetch gen2 cards: %s"' % (e))
            hdr = pyfits.Header()

        try:
            detectorTemp = self.actor.models['meb'].keyVarDict['temps'].valueList[1]
        except Exception as e:
            cmd.warn('text="FAILED to fetch MEB cards: %s"' % (e))
            detectorTemp = str(np.nan)

        try:
            config = pfsConfig.getConfigDict('mcs')
            detectorId = config['serial']
            gain = config['gain']
        except Exception as e:
            cmd.warn('text="FAILED to fetch config cards: %s"' % (e))
            detectorId = 'MCS cam#1'
            gain = 2.24
            
        hdr.append(('DETECTOR', detectorId, 'Name of the detector/CCD'))
        hdr.append(('GAIN', gain, '[e-/ADU] AD conversion factor'))
        hdr.append(('DET-TMP', detectorTemp, '[degC] Detector temperature'))
        try:
            instCards = self._getInstHeader(cmd)
            hdr.add_comment('Subaru Device Dependent Header Block for PFS-MCS')
            hdr.extend(instCards, bottom=True)
        except Exception as e:
            cmd.warn(f'text="FAILED to gather instrument cards: {e}"')

        return hdr

    def dumpCentroidtoDB(self, cmd):
        """Connect to database and return json string to an attribute."""
        
        conn = self.conn
        cmd.diag(f'text="dumping centroids to db {conn}"')
        
        # The section should be removed, replaced by the createTables command.
        if self.newTable:
            self._makeTables(conn, doDrop=True)
            cmd.inform('text="Centroid table created. "')
        else:
            #self._makeTables(conn, doDrop=False)
            cmd.inform('text="Attaching centroid to exsiting table. "')
        
        frameID=self.actor.exposureID
        buf = self._writeCentroids(self.centroids,1,frameID,1,conn)

        cmd.inform('text="Centroids of exposure ID %06d00 dummped. "'%(frameID))

    def _makeImageHeader(self, cmd):
        """ Create a complete WCS header.

        Notes
        ----
        - We need to get image center from config file.
        - So far, just spit out MCS-to-PFI linear conversion. Later, add SIP terms, and MCS-to-sky.
        - Needs to be pulled out into actorcore or pfs_utils.

        Returns
        -------
        hdr : pyfits.Header instance.

        """

        gen2Keys = self.actor.models['gen2'].keyVarDict

        imageCenter = (8960//2, 5778//2)
        rot = gen2Keys['tel_rot'][1]
        alt = gen2Keys['tel_axes'][1]
        az = gen2Keys['tel_axes'][0]

        cmd.debug(f'text="center={imageCenter} axes={az},{alt},{rot}"')

        hdr = None
        wcs, sip = pfsWcs.WCSParameters('mcs_pfi', imageCenter, rot, alt, az)
        hdr = wcs.to_header()

        return hdr

    def _doExpose(self, cmd, expTime, expType):
        """ Take an exposure and save it to disk. """

        nameQ = self.requestNextFilename(cmd)
        cmd.diag(f'text="new exposure"')
        if self.simulationPath is None:
            filename = 'scratchFile'
            image = self.actor.camera.expose(cmd, expTime, expType, filename, doCopy=False)
        else:
            image = self.getNextSimulationImage(cmd)
        cmd.diag(f'text="done: {image.shape}"')

        cmd.diag('text="reading filename"')
        try:
            filename = nameQ.get(timeout=5.0)
        except queue.Empty:
            cmd.warn('text="failed to get a new filename in time"')
            raise
        cmd.diag(f'text="read filename: {filename}"')

        hdr = self._constructHeader(cmd, filename, expType, expTime)
        cmd.diag(f'text="hdr done: {len(hdr)}"')
        phdu = pyfits.PrimaryHDU(header=hdr)

        try:
            imgHdr = self._makeImageHeader(cmd)
        except Exception as e:
            cmd.warn(f'text="FAILED to generate WCS header: {e}"')
            imgHdr = pyfits.Header()

        imgHdu = pyfits.CompImageHDU(image, name='IMAGE', compression_type='RICE_1')
        imgHdu.header.extend(imgHdr)
        hduList = pyfits.HDUList([phdu, imgHdu])

        # Patch core FITS card comments to match Subaru requirements.
        imgHdr = imgHdu.header
        imgHdr.set('BZERO', comment='Real=fits-value*BSCALE+BZERO')
        imgHdr.set('BSCALE', comment='Real=fits-value*BSCALE+BZERO')
        imgHdr.append(('BUNIT', 'ADU', 'Unit of original pixel values'))
        imgHdr.append(('BLANK', 32767, 'Value used for NULL pixels'))  # 12-bit camera
        imgHdr.append(('BIN-FCT1', 1, '[pixel] Binning factor of X axis'))
        imgHdr.append(('BIN-FCT2', 1, '[pixel] Binning factor of X axis'))

        pHdr =  phdu.header
        pHdr.set('EXTEND', comment='Presence of FITS Extension')

        hduList.writeto(filename, checksum=False, overwrite=True)
        cmd.inform('filename="%s"' % (filename))

        return filename, image
           
    def expose(self, cmd):
        """ Take an exposure. Optionally centroids. """

        cmdKeys = cmd.cmd.keywords
        doCentroid = 'doCentroid' in cmdKeys
        expType = cmdKeys[0].name
        if 'visit' in cmdKeys:
            self.actor.exposureID = cmdKeys['visit'].values[0]
        else:
            self.actor.exposureID = None

        if expType in ('bias', 'test'):
            expTime = self.expTime
        else:
            expTime = cmd.cmd.keywords["expTime"].values[0] * 1000

        if (expTime != self.expTime):
            self.actor.camera.setExposureTime(cmd,expTime)
 
        cmd.diag('text="Exposure time now is %d ms." '% (expTime))    
 
        filename, image = self._doExpose(cmd, expTime, expType)
        self.actor.image = image

        if doCentroid:
            self.runCentroid(cmd, doFinish=False)

        cmd.finish('exposureState=done')


    def doCentroidCoarse(self, cmd):
        

        pass
        
    def _doCentroid(self,cmd,image,fittype):

        pass

         
    def _encodeArray(self, array):
        """ Temporarily wrap a binary array transfer encoding. """

        # Actually, we want dtype,naxis,axNlen,base64(array)
        return base64.b64encode(array.tostring())

    def imageStats(self, cmd):
        
        doFinish = True
        cmd.inform('text="image median = %d." '% (np.median(self.actor.image))) 
        cmd.inform('text="image mean = %d." '% (self.actor.image.mean())) 
        cmd.inform('text="image min = %d." '% (self.actor.image.min())) 
        cmd.inform('text="image max = %d." '% (self.actor.image.max()))
        cmd.inform('text="image std = %d." '% (self.actor.image.std()))

        if doFinish:
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

    def setCentroidParams(self, cmd, doFinish=True):

        """

        Set the parameters for centroiding; placeholder for later routine to read from configuration file. 
        For each parameter it will check the command for keywords, and if it doesn't work, will go to a default.

        For the threshold, it will first check the command, then calculate from the image, then go to a default
        that won't crash the system. 

        """

        
        try:
            self.fhwm = cmd.cmd.keywords["fwhm"].values[0]
        except:
            self.fwhm = 3

        try:
            self.boxsize = cmd.cmd.keywords["boxsize"].values[0]
        except:
            self.boxsize = 3

        try:
            self.thresh = cmd.cmd.keywords["thresh"].values[0]
        except:
            try:
                self.thresh = self.actor.image.mean()+20*self.actor.image.std()
            except:
                self.thresh=2000
            
        try:
            self.rl = cmd.cmd.keywords["rl"].values[0]
        except:
            self.rl = -2.5

        try:
            self.rh = cmd.cmd.keywords["rh"].values[0]
        except:
            self.rh = 1.3

        try:
            self.sl = cmd.cmd.keywords["sl"].values[0]
        except:
            self.sl = 0.05

        try:
            self.sh = cmd.cmd.keywords["sh"].values[0]
        except:
            self.sh = 0.5

        if doFinish:
            cmd.finish('parameters=set')

    def createTables(self, cmd):
        """ Create SQL tables. """

        cmdKeys = cmd.cmd.keywords
        doDrop = 'drop' in cmdKeys
        suffix = cmdKeys['suffix'].values[0] if 'suffix' in cmdKeys else ''

        tables = dict()
        tables['mcsEngTable'] = '''id SERIAL PRIMARY KEY,
                                   datatime timestamp,
                                   frameId integer,
                                   moveId smallint,
                                   fiberId smallint,
                                   centroidX real, centroidY real,
                                   fwhmX real, fwhmY real,
                                   bgValue real, peakValue real'''

        with self.conn as conn:
            with conn.cursor() as curs:
                for t in tables:
                    tableName = f'{t}{suffix}'
                    if doDrop:
                        curs.execute(f'drop table if exists {tableName}')
                    createCmd = f'create table {tableName} ({tables[t]});'
                curs.execute(createCmd)

        cmd.finish('text="created tables')

    def runCentroid(self, cmd, doFinish=True):
        """ Measure centroids on the last acquired image. """

        cmdKeys = cmd.cmd.keywords
        self.newTable = "newTable" in cmdKeys
            
        cmd.debug('text="newTable value = %s"' % (self.newTable))

        image = self.actor.image
        
        cmd.inform(f'state="measuring cached image: {image.shape}"')
        a = centroid.centroid_only(image.astype('<i4'),
                                   self.fwhm, self.thresh, self.boxsize,
                                   2, self.sl, self.sh, self.rl, self.rh, 0)
        centroids=np.frombuffer(a,dtype='<f8')
        npoint=len(centroids)//7
        centroids=np.reshape(centroids,(npoint,7))
        centroids=centroids[:,0:6]

        self.centroids=centroids
            
        cmd.inform('text="%d centroids"'% (len(centroids)))
        cmd.inform('state="centroids measured"')
                        
        self.dumpCentroidtoDB(cmd)

        if doFinish:
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
        """ Create MCS tables. Someone said 20 measured things. """
        
        if self.simulationPath is None:
       
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
                    curs.execute('drop table if exists mcsPerFiber')
                    pass
                curs.execute(cmd)
        else:
            
            cmd = '''create table mcsEngTable (
                id SERIAL PRIMARY KEY,
                datatime timestamp,
                frameId integer,
                moveId smallint,
                fiberId smallint,
                
                centroidX real, centroidY real,
                fwhmX real, fwhmY real,
                bgValue real, peakValue real
                );'''
                
                    
            with conn.cursor() as curs:
                if doDrop:
                    curs.execute('drop table if exists mcsEngTable')
                curs.execute(cmd)
  
        conn.commit()
        
    def _writeCentroids(self, centArr, nextRowId, frameId, moveId, conn=None):
        """ Write all measurements for a given (frameId, moveId) """

        now = datetime.datetime.now()
        now.strftime("%Y-%m-%d %H:%M:%S")
            
        # Save measurements to a CSV buffer
        measBuf = io.StringIO()
        np.savetxt(measBuf, centArr, delimiter=',', fmt='%0.6g')
        measBuf.seek(0,0)

        # Let the database handle the primary key
        with conn.cursor() as curs:
            curs.execute("select * FROM mcsEngTable where false")
            colnames = [desc[0] for desc in curs.description]
        realcolnames = colnames[1:]

        buf = io.StringIO()
        for l_i in range(len(centArr)):
            line = '%s,%d,%d,%d,%s' % (now.strftime("%Y-%m-%d %H:%M:%S"), 
                                       frameId, moveId, l_i, measBuf.readline())
            buf.write(line)
        buf.seek(0,0)
            
        with conn:
            with conn.cursor() as curs:
                curs.copy_from(buf,'mcsEngTable',',',
                               columns=realcolnames)

        buf.seek(0,0)
        
        return buf

    def _readCentroids(self, conn, frameId, moveId):
        """ Read all measurements for a given (frameId, moveId)"""
        
        if self.simulationPath is None:
            buf = io.StringIO()
        
            cmd = """copy (select * from mcsPerFiber where frameId={frameId} and moveId={moveId}) to stdout delimiter ',' """
            with conn.cursor() as curs:
                curs.copy_expert(cmd, buf)
            conn.commit()
            buf.seek(0,0)
        
            # Skip the frameId, etc. columns.
            arr = np.genfromtxt(buf, dtype='f4',
                                delimiter=',',usecols=range(4,24))
        else:
            buf = io.StringIO()

            cmd = f"""copy (select * from mcsEngTable where frameId={frameId} and moveId={moveId}) to stdout delimiter ',' """
            with conn.cursor() as curs:
                curs.copy_expert(cmd, buf)
            conn.commit()
            buf.seek(0,0)

            # Skip the frameId, etc. columns.
            arr = np.genfromtxt(buf, dtype='f4',
                                delimiter=',',usecols=range(4,8))
        
        return arr

    def _writeDBtable(self):
        conn = psycopg2.connect("dbname='fps' user='pfs' host='localhost' password='pfspass'")
        cur = conn.cursor()
        cur.execute("select * from information_schema.tables where table_name=%s", ('test',))
        print(bool(cur.rowcount))
        conn.close()

    def timeTestFull(self,cmd):

        fwhm=3.
        hmin=2500
        boxsize=9
        expTime=1000.
        expType="object"

        #prelude stuff
        #fake table
        arr = np.random.uniform(0,1, (2394, 20)).astype('f4')

        #connect to database
        t1=time.time()
        for i in range(5):
            filename, image = self._doExpose(cmd, expTime, expType)
            self.actor.image = image.astype('<i4')

            a=centroid_only(self.actor.image,fwhm,hmin,boxsize)
            
            centroids=np.frombuffer(a,dtype='<f8')
            numpoints=len(centroids)//7
            cmd.inform('text="np = %d." '% (numpoints))
            
            self.actor.centroids=np.reshape(centroids,(numpoints,7))
            conn = psycopg2.connect("dbname='fps' user='pfs' host='localhost' password='pfspass'")
            self._makeTables(conn,doDrop=True)
            
            #db insertion
            buf = self._writeCentroids(arr,1,100,1,conn)
            cmd.inform('text="size = %d." '% (numpoints))
        t2=time.time()
        cmd.inform('text="time = %f." '% ((t2-t1)/1.))

        hduList.writeto(filename, checksum=False, overwrite=True)

        cmd.inform('filename="%s"' % (filename))

        return filename, image


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
