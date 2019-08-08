from builtins import range
from centroid import centroid_only
import pyfits
import numpy as np
import pylab as py
import os

image=pyfits.getdata('/Users/karr/Commissioning/Firsts/PFSC00102100.fits').astype('<i4')

fwhmx=4.
fwhmy=4.
thresh1=2500
thresh2=1200

boxFind=7
boxCent=6

nmin=10
nmax=90
maxIt=20


for i in range(100):
    a=centroid_only(image.astype('<i4'),fwhmx,fwhmy,thresh1,thresh2,boxFind,boxCent,nmin,nmax,maxIt,0)

    centroids=np.frombuffer(a,dtype='<f8')
    centroids=np.reshape(centroids,(len(centroids)//7,7))

#py.plot(centroids[:,0],centroids[:,1],'dg')
#py.show()

#print(centroids[:,0])

#for i in range(centroids.shape[0]):
#    print(centroids[i,0],centroids[i,1])

print(centroids.shape)
