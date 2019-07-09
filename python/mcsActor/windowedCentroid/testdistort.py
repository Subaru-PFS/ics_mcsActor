import numpy as np
import pyfits as pf
import pylab as py
import centroid as centroid
from heapq import nsmallest
import cv2
from numpy.random import rand
from numpy.random import randn


x=rand(200)-0.5
y=rand(200)-0.5


theta=0.5
s=1.2

a=s*np.cos(theta)
b=-s*np.sin(theta)
c=0.3
d=-0.1


x1=a*x-b*y+c+randn(len(x))*0.02
y1=b*x+a*y+d+randn(len(x))*0.02

print(x1.shape)

nn=len(x1)
pts1=np.zeros((1,nn,2))
pts2=np.zeros((1,nn,2))
for i in range(nn):
    pts1[0,i,0]=x[i]
    pts1[0,i,1]=y[i]
    pts2[0,i,0]=x1[i]
    pts2[0,i,1]=y1[i]
    #pts1[0,i,2]=0
    #pts2[0,i,2]=0

trans=cv2.estimateRigidTransform(pts1,pts2,False)
print(trans)
sx=np.sqrt(trans[0,0]**2+trans[0,1]**2)
sy=np.sqrt(trans[1,0]**2+trans[1,1]**2)
print("sx=",sx," sy=",sy)

trans[0,0]/=sx
trans[0,1]/=sx
trans[1,0]/=sx
trans[1,1]/=sx

points=cv2.transform(pts1,trans)

py.clf()
#py.scatter(pts1[0,:,0],pts1[0,:,1])
#py.scatter(pts2[0,:,0],pts2[0,:,1],color='g')
#py.scatter(points[0,:,0],points[0,:,1],color='r')
py.show()

py.clf()
py.quiver(pts1[0,:,0],pts1[0,:,1],pts2[0,:,0]-points[0,:,0],pts2[0,:,1]-points[0,:,1])
#py.show()


print(trans)
