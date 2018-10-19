import vis_plotting as visplot        #plotting routines 
import vis_calculations as viscalc    #calculation routines
import vis_coordinates as viscoords   #routines related to coordinate transforms
import numpy as np
from importlib import reload  #for development/debugging purposes
import matplotlib.pylab as plt


file1="MCST_010_001.fits"
file2="MCST_010_002.fits"
#file2="MCST_010_001.fits"

fwhm=3        
boxsize=9
thresh=2200
rl=-2.5
rh=1.5
sl=0.05
sh=0.5

inter=0
prefix="testrot"

image1=visplot.getImage(file1)
image2=visplot.getImage(file2)

x1,y1,fx1,fy1,back1,peak1=viscalc.getCentroids(image1,fwhm,boxsize,thresh,rl,rh,sl,sh)
x2,y2,fx2,fy2,back2,peak2=viscalc.getCentroids(image2,fwhm,boxsize,thresh,rl,rh,sl,sh)

visplot.checkCentroids(x1,y1,0,prefix,inter)
visplot.checkCentroids(x2,y2,0,prefix,inter)

reload(visplot)
reload(viscalc)
xt=3000.
yt=4500.
theta1=0.0
theta2=0.6

#x1_rot,y1_rot=viscalc.rotatePoints(x1,y1,xt,yt,theta1)
#x2_rot,y2_rot=viscalc.rotatePoints(x2,y2,xt,yt,theta2)

x1_rot,y1_rot=viscoords.transformGeneral(x1,y1,0,0,theta1,1,xt,yt)
x2_rot,y2_rot=viscoords.transformGeneral(x2,y2,0,0,theta2,1,xt,yt)

plt.clf()
plt.plot(x1_rot,y1_rot,'dg')
plt.plot(x2_rot,y2_rot,'ob')
plt.show()

#visplot.checkCentroids(x1_rot,y1_rot,0,prefix,inter)
#visplot.checkCentroids(x2_rot,y2_rot,0,prefix,inter)
xt=3000.
yt=4500.
theta=0.6

transformation,xd,yd,sx,sy,rotation = viscalc.matchRot(x1_rot,y1_rot,x2_rot,y2_rot,xt,yt,theta)
print(xd,yd,sx,sy,rotation)
print(transformation)

xt,yt=viscalc.getCentre(transformation)
theta=rotation
transformation,xd,yd,sx,sy,rotation = viscalc.matchRot(x1_rot,y1_rot,x2_rot,y2_rot,xt,yt,theta)

xt,yt=viscalc.getCentre(transformation)


#visplot.checkMatched(xx,yy,x1_rot,y1_rot,prefix,inter)
