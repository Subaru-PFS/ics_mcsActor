from builtins import range
from centroid import get_homes_call
from centroid import centroid_coarse_call
from centroid import centroid_fine_call
from centroid import centroid_only
import pyfits
import numpy as np
import pylab as py

image=pyfits.getdata('./MCSA0000000010.fits').astype('<i4')

fwhm=3.
hmin=3000
boxsize=9

a=centroid_only(image, fwhm, hmin, boxsize)

centroids=np.frombuffer(a,dtype='<f8')
centroids=np.reshape(centroids,(len(centroids)//7,7))

print(centroids.shape)
#py.plot(centroids[:,0],centroids[:,1],'dg')
#py.show()

#print(centroids[:,0])

#for i in range(centroids.shape[0]):
#    print(centroids[i,0],centroids[i,1])

print(centroids.shape)
