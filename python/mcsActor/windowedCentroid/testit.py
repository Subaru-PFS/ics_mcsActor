from centroid import centroid_only
from astropy.io.fits import getdata
import numpy as np
import matplotlib.pyplot as plt
import os
image=getdata('/Users/karr/Science/PFS/Firsts/PFSC00366600.fits')


os.system('date')
for i in range(50):
    a=centroid_only(np.ascontiguousarray(image.astype('<i4')[2000:6000,2500:7000]),0,0,2200,1200,10,6,10,90,20,0)
    #a=centroid_only(image.astype('<i4'),0,0,2200,1200,10,6,10,90,20,0)
    centroids=np.frombuffer(a,dtype='<f8')
    centroids=np.reshape(centroids,(len(centroids)//7,7))
    
    #ind=np.where(centroids[:,1] < 1000)
    #print(centroids[ind,0],centroids[ind,1])

    #print(centroids[:,0].min(),centroids[:,0].max(),centroids[:,1].min(), centroids[:,1].max())
os.system('date')
print(centroids.shape)
#plt.ion()
#fig,ax=plt.subplots()
#ax.scatter(centroids[:,0],centroids[:,1])
#fig.show()
#plt.savefig("cent.png")

