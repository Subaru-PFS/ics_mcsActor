from builtins import range
from centroid import get_homes_call
from centroid import centroid_coarse_call
from centroid import centroid_fine_call
import astropy.io.fits as pyfits

import numpy as np
import pylab as py

#read in home image
image=pyfits.getdata('./home.fits').astype('<i4')

#get homes positions
a=get_homes_call(image)

#turn into numpy array
homes=np.frombuffer(a,dtype='<f8')

#get first positoin/arc

image=pyfits.getdata('./first_move.fits').astype('<i4')
arc_image=pyfits.getdata('./first_move_arc.fits').astype('<i4')

#call the routine
b=centroid_coarse_call(image,arc_image,homes)

#into numpy array
homepos=np.frombuffer(b,dtype=[('xp','<f8'),('yp','<f8'),('xt','<f8'),('yt','<f8'),('xc','<f8'),('yc','<f8'),('x','<f8'),('y','<f8'),('peak','<f8'),('back','<f8'),('fx','<f8'),('fy','<f8'),('qual','<f4'),('idnum','<f4')])

#again for the secon dmove
image=pyfits.getdata('./second_move.fits').astype('<i4')
arc_image=pyfits.getdata('./second_move.fits').astype('<i4')

b=centroid_coarse_call(image,arc_image,homes)

homepos=np.frombuffer(b,dtype=[('xp','<f8'),('yp','<f8'),('xt','<f8'),('yt','<f8'),('xc','<f8'),('yc','<f8'),('x','<f8'),('y','<f8'),('peak','<f8'),('back','<f8'),('fx','<f8'),('fy','<f8'),('qual','<f4'),('idnum','<f4')])

#now the third position
image=pyfits.getdata('./third_move.fits').astype('<i4')

npos=homepos.shape[0]

xp=np.zeros((npos))
yp=np.zeros((npos))

for i in range(npos):
	xp[i]=homepos[i][6]
	yp[i]=homepos[i][7]

c=centroid_fine_call(image,homes,xp,yp)

homepos=np.frombuffer(c,dtype=[('xp','<f8'),('yp','<f8'),('xt','<f8'),('yt','<f8'),('xc','<f8'),('yc','<f8'),('x','<f8'),('y','<f8'),('peak','<f8'),('back','<f8'),('fx','<f8'),('fy','<f8'),('qual','<f4'),('idnum','<f4')])


#uncomment this to do some plots to check results

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
#	xp[i]=homepos[i][0]
#	yp[i]=homepos[i][1]
#	xc[i]=homepos[i][4]
#	yc[i]=homepos[i][5]
#	xh[i]=homepos[i][6]
#	yh[i]=homepos[i][7]
#
#
#py.clf()
#py.plot(xp,xc,'dg')
#py.show()


#py.clf()
#py.plot(xc,yc,'dg')
#py.plot(xh,yh,'ob')
#fig=py.gcf()
#ax=fig.gca()
#
#
#for i in range(len(xc)):
#	circle=py.Circle((xh[i],yh[i]),55,color="black",fill=False)
#	ax.add_artist(circle)
#
#	py.plot([xh[i],xc[i]],[yh[i],yc[i]])
#
#
#py.show()
