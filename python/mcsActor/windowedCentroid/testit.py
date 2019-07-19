from centroid import centroid_only
from astropy.io.fits import getdata
import numpy as np

image=getdata('/Users/karr/Science/PFS/Firsts/PFSC00366600.fits')
a=centroid_only(image.astype('<i4'),0,0,2200,1200,10,6,10,90,20,1)
print(len(a))
centroids=np.frombuffer(a,dtype='<f8')
centroids=np.reshape(centroids,(len(centroids)//7,7))

import matplotlib.pyplot as plt

plt.ion()
fig,ax=plt.subplots()
ax.scatter(centroids[:,0],centroids[:,1])
fig.show()

ind=np.where(centroids[:,1] < 1000)
print(centroids[ind,0],centroids[ind,1])
