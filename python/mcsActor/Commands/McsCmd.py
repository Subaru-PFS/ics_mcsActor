#!/usr/bin/env python

import time
import datetime
import io
import pathlib
import queue
import threading
import sep

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

import mcsActor.windowedCentroid.centroid as centroid
import mcsActor.Visualization.mcsRoutines as mcsTools

import pandas as pd

import numpy as np

# matplotlib.use('Agg')
Bool=bool

class McsCmd(object):

    def __init__(self, actor):
        # This lets us access the rest of the actor.
        self.actor = actor
        self.expTime = 1000
        self.newTable = None
        self.simulationPath = None
        self._conn = None

        self.findThresh = None
        self.centThresh =  None

        self.db='db-ics'
        self.setCentroidParams(None)
        
        # Declare the commands we implement. When the actor is started
        # these are registered with the parser, which will call the
        # associated methods when matched. The callbacks will be
        # passed a single argument, the parsed and typed command.
        #
        self.vocab = [
            ('ping', '', self.ping),
            ('status', '', self.status),
            ('expose', '@(bias|test) [<frameId>]', self.expose),
            ('expose', '@(dark|flat) <expTime> [<frameId>]', self.expose),
            ('expose', '@object <expTime> [<frameId>] [@noCentroid] [@doCentroid] [@doFibreID] [@simDot]', self.expose),
            ('runCentroid', '[@newTable]', self.runCentroid),
            #('runFibreID', '[@newTable]', self.runFibreID),
            ('reconnect', '', self.reconnect),
            ('resetThreshold','',self.resetThreshold),
            ('setCentroidParams','[<fwhmx>] [<fwhmy>] [<boxFind>] [<boxCent>] [<nmin>] [<nmax>] [<maxIt>]',
             self.setCentroidParams),
            ('calcThresh','[<threshMethod>] [<threshSigma>] [<threshFact>]', self.calcThresh),
            ('simulate', '<path>', self.simulateOn),
            ('simulate', 'off', self.simulateOff),
        ]

        # Define typed command arguments for the above commands.
        self.keys = keys.KeysDictionary("mcs_mcs", (1, 1),
                                        keys.Key("expTime", types.Float(), help="The exposure time, seconds"),
                                        keys.Key("expType", types.String(), help="The exposure type"),
                                        keys.Key("frameId", types.Int(), help="exposure frameID"),
                                        keys.Key("path", types.String(), help="Simulated image directory"),
                                        keys.Key("getArc", types.Int(), help="flag for arc image"),
                                        keys.Key("fwhmx", types.Float(), help="X fwhm for centroid routine"),
                                        keys.Key("fwhmy", types.Float(), help="Y fwhm for centroid routine"),
                                        keys.Key("boxFind", types.Int(), help="box size for finding spots"),
                                        keys.Key("boxCent", types.Int(), help="box size for centroiding spots"),
                                        keys.Key("nmin", types.Int(), help="minimum number of points for spot"),
                                        keys.Key("nmax", types.Int(), help="max number of points for spot"),
                                        keys.Key("maxIt", types.Int(), help="maximum number of iterations for centroiding"),
                                        keys.Key("findSigma", types.Float(), help="threshhold for finding spots"),
                                        keys.Key("centSigma", types.Float(), help="threshhold for calculating moments of spots"),
                                        keys.Key("matchRad", types.Int(), help="radius in pixels for matching positions"),
                                        keys.Key("threshMethod", types.String(), help="method for thresholding"),
                                        keys.Key("threshSigma", types.Float(), help="simga for sigma-clipped RMS of image"),
                                        keys.Key("threshFact", types.Float(), help="factor for thresholding"),
                                        keys.Key("fieldID", types.String(), help="fieldID for getting instrument parameters")
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
            connString = "dbname='opdb_asrd' user='pfs' host="+self.db+" password="+passstring
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
        self.actor.connectCamera(cmd, camera='rmod')
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

    def requestNextFilename(self, cmd, frameId):
        """ Return a queue which will eventually contain a filename. """

        q = queue.Queue()

        def worker(q=q, cmd=cmd):
            filename = self.getNextFilename(cmd, frameId)
            q.put(filename)
            q.task_done()

        task = threading.Thread(target=worker,
                                name='visitFetcher', daemon=True)
        task.start()

        return q


    def _insertPFSVisitID(self, visitid):
        with self.conn:
            with self.conn.cursor() as curs:
                postgres_insert_query = """ INSERT INTO pfs_visit (pfs_visit_id, pfs_visit_description) VALUES (%s,%s)"""
                record_to_insert = (visitid, "MCS exposure")
                curs.execute(postgres_insert_query, record_to_insert)


    def getNextFilename(self, cmd, frameId):
        """ Fetch next image filename. 

        In real life, we will instantiate a Subaru-compliant image pathname generating object.  

        """

        if frameId is None:
            ret = self.actor.cmdr.call(actor='gen2', cmdStr='getVisit', timeLim=10.0)
            if ret.didFail:
                raise RuntimeError("getNextFilename failed getting a visit number in 10s!")

            visit = self.actor.models['gen2'].keyVarDict['visit'].valueList[0]
            frameId = visit * 100
            self._insertPFSVisitID(visit)

        path = os.path.join("$ICS_MHS_DATA_ROOT", 'mcs', time.strftime('%Y-%m-%d', time.gmtime()))
        path = os.path.expandvars(os.path.expanduser(path))
        if not os.path.isdir(path):
            os.makedirs(path, 0o755)

        newpath = os.path.join(path, 'PFSC%08d.fits' % (frameId))

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

    def dumpCentroidtoDB(self, cmd, frameId):
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
        
        buf = self._writeCentroids(self.centroids,1,frameId,1,conn)

        cmd.inform('text="Centroids of exposure ID %08d dumped."' % (frameId))

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

    def _doExpose(self, cmd, expTime, expType, frameId, mask=None):
        """ Take an exposure and save it to disk. """

        nameQ = self.requestNextFilename(cmd, frameId)
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
        
        if mask is not None:
            cmd.inform(f'text="mask image shape: {mask.shape} type:{mask.dtype}"')
            cmd.inform(f'text="image shape: {image.shape} type:{image.dtype}"')
            image = image*mask.astype('uint16')
            
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
    
    def resetThreshold(self,cmd):
        self.findThresh = None
        self.centThresh =  None
        cmd.finish('Centroid threshold=d=reset')
       
    def expose(self, cmd):
        """ Take an exposure. Optionally centroids. Optionally FibreID """

        cmdKeys = cmd.cmd.keywords

        # Switch from default no centroids to default do centroids
        noCentroidArg = 'noCentroid' in cmdKeys
        if noCentroidArg:
            doCentroid = False
        else:
            doCentroid = True
        
        doFibreID = 'doFibreID' in cmdKeys

        simDotArg = 'simDot' in cmdKeys
        if simDotArg:
            simDot = True
            dotmask = self._makeDotMask(cmd)
        else:
            simDot = False



        cmd.inform('text="doCentroid = %s." '%{doCentroid})
        expType = cmdKeys[0].name
        if 'frameId' in cmdKeys:
            frameId = cmdKeys['frameId'].values[0]
        else:
            frameId = None

        if expType in ('bias', 'test'):
            expTime = self.expTime
        else:
            expTime = cmd.cmd.keywords["expTime"].values[0] * 1000

        if (expTime != self.expTime):
            self.actor.camera.setExposureTime(cmd,expTime)
        
        cmd.inform('text="Exposure time now is %d ms." '% (expTime))    
        if simDot is True:
            filename, image = self._doExpose(cmd, expTime, expType, frameId, mask=dotmask)
        else:
            filename, image = self._doExpose(cmd, expTime, expType, frameId)

        if frameId is None:
            filename = pathlib.Path(filename)
            frameId = int(filename.stem[4:], base=10)
        self.actor.image = image

        self.handleTelescopeGeometry(cmd, filename, frameId, expTime)

        if doCentroid:
            cmd.inform('text="Setting centroid parameters." ')
            self.setCentroidParams(cmd)
   
            #self.calcThresh(cmd)
            if self.findThresh is None:
                cmd.inform('text="Calculating threshold." ')
                self.calcThresh(cmd)

            cmd.inform('text="Running centroid on current image" ')
            self.runCentroidSEP(cmd)
            
            #if (self.nCentroid < 2000):
            #    self.calcThresh(cmd)
            #    self.runCentroid(cmd)       
            
            cmd.inform('text="Sending centroid data to database" ')
            self.dumpCentroidtoDB(cmd, frameId)
        
        if doFibreID:
            self.runFibreID(cmd, doFinish=False)
            self.dumpCentroidtoDB(cmd, frameId)


        cmd.finish('exposureState=done')

    def handleTelescopeGeometry(self, cmd, filename, frameId, expTime):

        if self.simulationPath is None:
            # We are live: use Gen2 telescope info.
            gen2Model = self.actor.models['gen2'].keyVarDict

            axes = gen2Model['tel_axes'].getValue()
            az, alt = axes

            rot = gen2Model['tel_rot'].getValue()
            posAngle, instrot = rot
            now = datetime.datetime.now()
        else:
            # We are reading images from disk: get the geometry from the headers.
            hdr = fitsio.read_header(str(filename), 0)
            simPath = hdr['W_MCSMNM']
            simHdr = fitsio.read_header(str(simPath), 0)
            cmd.inform('text="loaded telescope info from %s"'% (simPath))
            az = simHdr['AZIMUTH']
            alt = simHdr['ALTITUDE']
            expTime = simHdr['EXPTIME']
            instrot = simHdr['INR-STR']

        now = datetime.datetime.now()

        cmd.inform(f'text="az={az} alt={alt} instrot={instrot}"')
        # Packing information into data structure
        telescopeInfo = {'frameid': frameId,
                         'visitid': self.actor.models['gen2'].keyVarDict['visit'].valueList[0],
                         'starttime': now.strftime("%Y-%m-%d %H:%M:%S"),
                         'exptime': expTime,
                         'altitude': alt,
                         'azimuth': az,
                         'instrot': instrot}
        print("telescopeInfo")
        self._writeTelescopeInfo(cmd,telescopeInfo, self.conn)

 

    def calcThresh(self, cmd):

        """  Calculate thresholds for finding/centroiding from image 

        3 methods: fieldID uses the known system geometry to figure out the right region
                   calib is for calibration when the system is not known, calculates from the image characteristics
                   direct sets teh values manually (backup method)


        """
        image = self.actor.image

        self.findThresh,self.centThresh,xrange,yrange = mcsTools.getThresh(image,
            'calib',4,2,self.findSigma,self.centSigma)
        
        cmd.inform('text="findThresh = %d, centThresh = %d." '%(self.findThresh,self.centThresh))    
 



        
    def setCentroidParams(self, cmd):

        """

        Set the parameters for centroiding; placeholder for later routine to read from configuration file. 
        For each parameter it will check the command for keywords, and if it doesn't work, will go to a default.

        For the threshold, it will first check the command, then calculate from the image, then go to a default
        that won't crash the system. 

        """

        
        try:
            self.fhwmx = cmd.cmd.keywords["fwhmx"].values[0]
        except:
            self.fwhmx = 0

        try:
            self.fhwmy = cmd.cmd.keywords["fwhmy"].values[0]
        except:
            self.fwhmy = 0

        try:
            self.boxFind = cmd.cmd.keywords["boxFind"].values[0]
        except:
            self.boxFind = 10
            
        try:
            self.boxCent = cmd.cmd.keywords["boxCent"].values[0]
        except:
            self.boxCent = 6

        try:
            self.findSigma = cmd.cmd.keywords["findSigma"].values[0]
        except:
            self.findSigma = 35
            
        try:
            self.centSigma = cmd.cmd.keywords["centSigma"].values[0]
        except:
            self.centSigma = 15

        try:
            self.nmin = cmd.cmd.keywords["nmin"].values[0]
        except:
            self.nmin = 10
            
        try:
            self.nmax = cmd.cmd.keywords["nmax"].values[0]
        except:
            self.nmax = 90
            
        try:
            self.maxIt = cmd.cmd.keywords["maxIt"].values[0]
        except:
            self.maxIt = 20
            
        try:
            self.matchRad = cmd.cmd.keywords["matchRad"].values[0]
        except:
            self.matchRad = 20

    def runCentroidSEP(self, cmd):
        cmdKeys = cmd.cmd.keywords
        self.newTable = "newTable" in cmdKeys
            
        cmd.debug('text="newTable value = %s"' % (self.newTable))

        image = self.actor.image
        
        cmd.inform(f'state="measuring cached image: {image.shape}"')
        centroids = sep.extract(image.astype(float), 300)
        npoint=centroids.shape[0]
        tCentroids = np.zeros((npoint,8))

        tCentroids[:,0]=np.arange(npoint)+1
        tCentroids[:,1]=centroids['x']
        tCentroids[:,2]=centroids['y']
        tCentroids[:,3]=centroids['x2']
        tCentroids[:,4]=centroids['y2']
        tCentroids[:,5]=centroids['xy']
        tCentroids[:,6]=centroids['thresh']
        tCentroids[:,7]=centroids['peak']

        self.centroids=tCentroids
        self.nCentroid = len(centroids)    
        cmd.inform('text="%d centroids"'% (len(centroids)))
        cmd.inform('state="centroids measured"')

    def runCentroid(self, cmd):
        """ Measure centroids on the last acquired image. """

        cmdKeys = cmd.cmd.keywords
        self.newTable = "newTable" in cmdKeys
            
        cmd.debug('text="newTable value = %s"' % (self.newTable))

        image = self.actor.image
        
        cmd.inform(f'state="measuring cached image: {image.shape}"')
        a = centroid.centroid_only(image.astype('<i4'),
                                   self.fwhmx, self.fwhmy, self.findThresh, self.centThresh, self.boxFind, self.boxCent, 
                                   self.nmin, self.nmax, self.maxIt, 0)
        centroids=np.frombuffer(a,dtype='<f8')
        centroids=np.reshape(centroids,(len(centroids)//7,7))
    
        npoint=centroids.shape[0]
    
    
        tCentroids = np.zeros((npoint,8))
        tCentroids[:,1:]=centroids
    
        bg = np.copy(tCentroids[:,6])
        peak = np.copy(tCentroids[:,5])

        tCentroids[:,6] = peak
        tCentroids[:,5] = bg
        
        self.centroids=tCentroids
        self.nCentroid = len(centroids)    
        cmd.inform('text="%d centroids"'% (len(centroids)))
        cmd.inform('state="centroids measured"')
                        
    
    
    def _writeTelescopeInfo(self, cmd, telescopeInfo, conn = None):

        # Let the database handle the primary key
        with conn:
            with conn.cursor() as curs:
                curs.execute('select * FROM "mcs_exposure" where false')
                colnames = [desc[0] for desc in curs.description]
            realcolnames = colnames[0:]
        
        colname = []
        for i in realcolnames:
            x='"'+i+'"'
            colname.append(x)
        
        buf = io.StringIO()

        """
          TODO: Those are the fake values for making PFI to work now, adding actual code later
        """
        adc_pa = 0
        dome_temperature = 5
        dome_pressure = 101
        dome_humidity = 0.3
        outside_temperature = 5
        outside_pressure = 101
        outside_humidity = 0.3        
        mcs_cover_temperature = 5
        mcs_m1_temperature =6
        taken_at = telescopeInfo['starttime']
        taken_in_hst_at = telescopeInfo['starttime']
        line = '%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%s,%s' % (telescopeInfo['frameid'], 
            telescopeInfo['visitid'], telescopeInfo['exptime']/1000.0,
            telescopeInfo['altitude'],telescopeInfo['azimuth'],telescopeInfo['instrot'],
            adc_pa,dome_temperature,dome_pressure,dome_humidity,outside_temperature,outside_pressure,
            outside_humidity,mcs_cover_temperature,mcs_m1_temperature,taken_at,taken_in_hst_at)

        buf.write(line)
        buf.seek(0,0)
            
        with conn:
            with conn.cursor() as curs:
                curs.copy_from(buf,'"mcs_exposure"',',',
                               columns=colname)

        buf.seek(0,0)
        
        cmd.inform('text="Telescope information for frame %s populated."' % (telescopeInfo['frameid']))

        return buf


    def _writeCentroids(self, centArr, nextRowId, frameId, moveId, conn=None):
        """ Write all measurements for a given (frameId, moveId) """

        now = datetime.datetime.now()
        now.strftime("%Y-%m-%d %H:%M:%S")
            
        # Save measurements to a CSV buffer
        measBuf = io.StringIO()
        

        data = np.insert(centArr, 5, 0, axis=1)
       
        np.savetxt(measBuf, data[:,1:8], delimiter=',', fmt='%0.6g')
        measBuf.seek(0,0)

        # Let the database handle the primary key
        with conn:
            with conn.cursor() as curs:
                curs.execute('select * FROM "mcs_data" where false')
                colnames = [desc[0] for desc in curs.description]
            realcolnames = colnames[:]
        
        colname = []
        for i in realcolnames:
            x='"'+i+'"'
            colname.append(x)
        buf = io.StringIO()
        for l_i in range(len(centArr)):
            #line = '%s,%d,%d,%d,%s' % (now.strftime("%Y-%m-%d %H:%M:%S"), 
            #                           frameId, moveId, l_i+1, measBuf.readline())
            line = '%d,%d,%s' % (frameId, l_i+1, measBuf.readline())
            buf.write(line)
        buf.seek(0,0)
            
        with conn:
            with conn.cursor() as curs:
                curs.copy_from(buf,'"mcs_data"',',',
                               columns=colname)

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

            cmd = f"""copy (select * from 'mcsData' where frameId={frameId} and moveId={moveId}) to stdout delimiter ',' """
            with conn.cursor() as curs:
                curs.copy_expert(cmd, buf)
            conn.commit()
            buf.seek(0,0)

            # Skip the frameId, etc. columns.
            arr = np.genfromtxt(buf, dtype='f4',
                                delimiter=',',usecols=range(4,8))
        
        return arr

   
    def _makeDotMask(self,cmd):
        """ Make dot mask for simulation or proessing purpose)"""

        # read dot location
        dotfile = '/home/pfs/mhs/devel/pfs_instdata/data/pfi/dot/dot_measurements_20210428_el90_rot+00_ave.csv'    
        dotpos = pd.read_csv(dotfile)
        

        unit_height = 30
        unit_weight = 30
        r_mean=np.around(np.mean(dotpos['r_tran'].values)).astype('int')
        cmd.inform('text="Making dot image"')

        mask = self._create_circular_mask(unit_height, unit_weight, radius=10)

        masked_img = np.zeros([unit_height,unit_weight])+1

        masked_img[mask] = 0
        # To-do: change here for the image size
        dotmask = np.zeros([7096,10000])+1

        for i in range(len(dotpos)):
            xstart=np.around(dotpos['x_tran'].values[i]-(unit_weight/2)).astype('int')

            ystart=np.around(dotpos['y_tran'].values[i]-(unit_height/2)).astype('int') 

            dotmask[ystart:ystart+unit_height,xstart:xstart+unit_weight] = masked_img

        return dotmask
    
    def _create_circular_mask(self, h, w, center=None, radius=None):

        if center is None: # use the middle of the image
            center = (int(w/2), int(h/2))
        if radius is None: # use the smallest distance between the center and image walls
            radius = min(center[0], center[1], w-center[0], h-center[1])

        Y, X = np.ogrid[:h, :w]
        dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

        mask = dist_from_center <= radius
        return mask        