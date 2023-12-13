
import numpy as np
from multiprocessing import Queue
from multiprocessing import Process
from multiprocessing import shared_memory

import astropy.io.fits as pyfits
import sep
import time
import logging
from os import environ

class speedCentroid(object):
    def __init__(self, image):
        logging.basicConfig(format="%(asctime)s.%(msecs)03d %(levelno)s %(name)-10s %(message)s",
                            datefmt="%Y-%m-%dT%H:%M:%S")
        self.logger = logging.getLogger('speedCentroid')
        self.logger.setLevel(logging.INFO)
    
        self.image = image
        self.cores = 12
        
        self.logger.info(f'Estimating background')
        bkg = sep.Background(image[3200:3400,5200:5400].astype(float), bw=64, bh=64, fw=3, fh=3)
        
        self.logger.info(f'Setting background to be {bkg.globalrms}')
        self.bkgglobalrms = bkg.globalrms

        

        self.logger.info(f'Creating share memory')
        self.shm_name, self.shm_dtype, self.shm_shape = self.shareData(image)
        
        
        self.logger.info(f"{self.shm_name} {self.shm_dtype} {self.shm_shape}")
        
    def shareData(self,data):
        
        shm = shared_memory.SharedMemory(name=None, create=True, size=data.nbytes)
        dataBuf = np.ndarray(dtype=data.dtype, shape=data.shape, buffer=shm.buf)
        dataBuf[:] = data[:]
        shm.close()
        return (shm.name, dataBuf.dtype.name, dataBuf.shape)
    
    def runCentroid(self, image, queue=None, yshift=None):

        centroids = sep.extract(image.astype(float), 20 , err=self.bkgglobalrms,
            filter_type='conv', minarea=10)
        
        if yshift is not None:
            centroids['y'] += yshift 
        
        if queue is not None:
            queue.put(centroids)
    
    def arrangeCentroid(self):
        all_sources = self.q[0]

        for i in range(1, self.cores):            
            all_sources = np.concatenate((all_sources, self.q[i]), axis=0)

        self.centroids = all_sources
    
    
    def runCentroidMP(self):
        shm = shared_memory.SharedMemory(name=self.shm_name, create=False)
        _data = np.ndarray(dtype=self.shm_dtype, shape=self.shm_shape, buffer=shm.buf)
        
        yshift=self.image.shape[0]// self.cores
        q = Queue()
        
        process = list()
        
        for i in range(self.cores):
            self.logger.info(f'image size = {_data[yshift*i:yshift*(i+1),:].shape}')
            x = Process(target=self.runCentroid, args=(_data[yshift*i:yshift*(i+1),:],),kwargs={'queue': q,'yshift': yshift*i})
            x.start()
            process.append(x)
        
        self.process = process
        
        self.q = [q.get() for _ in process] 
        for p in process:
            p.join()
    
    def close(self):
        shm = shared_memory.SharedMemory(name=self.shm_name, create=False)

        shm.close()
        shm.unlink()
        del shm


def main(argv=None):
    for i in range(10):

        file='/data/raw/2023-10-03/mcs/PFSC10055606.fits'
        image = pyfits.open(file)[1].data
        
        t0 = time.time()

        spCenMT = speedCentroid(image)
        spCenMT.runCentroidMP()
        spCenMT.arrangeCentroid()
        t1 = time.time()

        print(f'Multi-Process = {t1 - t0}')

        bkg = sep.Background(image.astype(float), bw=64, bh=64, fw=3, fh=3)
        centroids = sep.extract(image.astype(float), 20 , err=bkg.globalrms,
                filter_type='conv', minarea=10)

        t2 = time.time()

        print(f'Single-Process = {t2 - t1}')
        spCenMT.close()
    
if __name__ == "__main__":
    main()