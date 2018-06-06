

import numpy as np
import pyfits as pf
import pylab as py
import centroid as centroid
from heapq import nsmallest
import cv2
from matplotlib.cm import plasma
from scipy.stats import sigmaclip


def getImage(filename):

    image=pf.getdata(filename)
    return image

def getCentroids(image,fwhm,boxsize,thresh,scale,rl,rh,sl,sh):

    #centroid 
    a=centroid.centroid_only(image.astype('<i4'),fwhm,thresh,boxsize,2,sl,sh,rl,rh,0)
    centroids=np.frombuffer(a,dtype='<f8')
    npoint=len(centroids)//7
    centroids=np.reshape(centroids,(npoint,7))
    
    #extract useful values

    x=centroids[:,0]
    y=centroids[:,1]
    
    fx=(centroids[:,2])
    fy=(centroids[:,3])
    back=(centroids[:,4])
    peak=(centroids[:,5])

    return x,y,fx,fy,back,peak

def approximateCoords(x,y,region):

    ind=np.where((x > region[0]) & (x < region[1]) & (y > region[2]) & (y < region[3]))

    dd=x*x+y*y

    ind1=np.where(dd == dd.max())
    ind2=np.where(dd == dd.min())

    
    angle=np.arctan((y[ind1]-y[ind2])/(x[ind1]-x[ind2]))-np.pi/4.

    return x[ind2],y[ind2],angle

def maskinMM(smallBox):
                        
    #manually crop the close points
    xx=[]
    yy=[]

    #centre spots
    pp1=np.array([6,13,8,15])+0.5
    pp3=np.array([28,35,8,15])+0.5
    pp2=np.array([17,24,30,37])+0.5

    for i in range(43):
        for j in range(43):
            iin=1
            ii=i
            jj=j
            if((ii > pp1[0]) & (ii < pp1[1]) & (jj > pp1[2]) & (jj < pp1[3])):
               iin=0
            if((ii > pp2[0]) & (ii < pp2[1]) & (jj > pp2[2]) & (jj < pp2[3])):
               iin=0
            if((ii > pp3[0]) & (ii < pp3[1]) & (jj > pp3[2]) & (jj < pp3[3])):

               iin=0
            if((i==21) and (j==21)):
                iin=0
            if(iin==1):
               xx.append(ii*8)
               yy.append(jj*8)

    if(smallBox==1):
        for i in range(25):
            for j in range(25):
                x1=i*2+(pp1[0]+0.5)
                y1=j*2*(pp1[3]+0.5)
                x2=i*2+(pp2[0]+0.5)
                y2=j*2*(pp2[3]+0.5)
                x3=i*2+(pp3[0]+0.5)
                y3=j*2*(pp3[3]+0.5)

                xx.append(x1)
                xx.append(x2)
                xx.append(x3)

                yy.append(y1)
                yy.append(y2)
                yy.append(y3)

    return np.array(xx),np.array(yy)

def scaleCentroids(x,y,x1,y1,scale):

    xc=(x-x1)/scale
    yc=(y-y1)/scale

    return xc,yc

def scaleMask(xx,yy,angle,flip):
    #apply any rotation (for matching purposes only)
    xx=xx-168
    yy=yy-168

    if(flip==1):
        angle=angle+np.pi/2
        
    xnew=xx*np.cos(angle)-yy*np.sin(angle)
    ynew=xx*np.sin(angle)+yy*np.cos(angle)

        
    
    xx=xnew+168
    yy=ynew+168

    return xx,yy

def checkCentroids(xc,yc,xx,yy,cutrange,prefix):

    #quick plot of centroids
    py.clf()
    py.scatter(xc,yc)
    py.scatter(xx,yy,s=2,color='r')
    py.savefig("/Users/karr/"+prefix+"_checkpoints.jpg")

    if(cutrange==1):
        np.xlim([-8,344])
        np.ylim([-8,344])
        
def matchPoints(xc,yc,xx,yy,fx,fy,peak):    
    
    xs=np.zeros((len(xx)))
    ys=np.zeros((len(yy)))
    fxs=np.zeros((len(xx)))
    fys=np.zeros((len(xx)))
    peaks=np.zeros((len(xx)))
    
    #match points between two sets (nearest neighbour)
    #for i in range(len(xx)):
    for i in range(len(xx)):
        for j in range(len(xc)):

            rr=np.sqrt((xx[i]-xc[j])**2+(yy[i]-yc[j])**2)
            if(rr < 5):
                xs[i]=xc[j]
                ys[i]=yc[j]
                fxs[i]=fx[j]
                fys[i]=fy[j]
                peaks[i]=peak[j]

    return xs,ys,fxs,fys,peaks


def checkMatched(xx,yy,xs,ys,prefix):
    
    #quick plot of matched points
    py.clf()
    py.scatter(xs,ys)
    py.scatter(xx,yy,s=2,color='r')
    py.savefig("/Users/karr/"+prefix+"_checkpoints1.jpg")


def simpleDistortion(xx,yy,xs,ys):    

    #subtract average to align

    nn=len(xs)
    pts1=np.zeros((1,nn,2))
    pts2=np.zeros((1,nn,2))
    xav1=xs.mean()
    yav1=ys.mean()
    xav2=xx.mean()
    yav2=yy.mean()

    for i in range(len(xs)):
        pts1[0,i,0]=xs[i]-xav1
        pts1[0,i,1]=ys[i]-yav1
        pts2[0,i,0]=xx[i]-xav2
        pts2[0,i,1]=yy[i]-yav2

    #cv2 needs float32
    
    pts1=np.float32(pts1)
    pts2=np.float32(pts2)

    #get the affine transformation - translation, rotation, scaling
    #uses opencv library
    
    transformation = cv2.estimateRigidTransform(pts1, pts2, False)
    print(transformation)
 
    #calculate the scaling from the result
    sx=np.sqrt(transformation[0,0]**2+transformation[0,1]**2)
    sy=np.sqrt(transformation[1,0]**2+transformation[1,1]**2)
    
    print("sx=",sx," sy=",sy)

    #remove the scaling from the transformation matrix
    transformation[0,0]/=sx
    transformation[0,1]/=sx
    transformation[1,0]/=sy
    transformation[1,1]/=sy

    print(transformation)

    #calculate the transformed points
    points=cv2.transform(pts1,transformation)

    #difference between expected and actual point positions
    diffx=(pts2[0,:,0]-points[0,:,0])
    diffy=(pts2[0,:,1]-points[0,:,1])

    #in mm and in % of field
    
    c=np.sqrt(diffx*diffx+diffy*diffy)
    c1=np.sqrt(diffx*diffx+diffy*diffy)/336.*100

    return c,c1,pts1,pts2,diffx,diffy,

def plotDistortion(c,c1,pts1,pts2,diffx,diffy,fxs,fys,peaks,limit,prefix):

    if(limit > 0):
        #a quick cludge to filter out bad quality points in poorly focussed images
    
        ind=np.where((c < limit) & (fxs < 10) & (fys < 10) & (fxs > 0) & (fys > 0) & (peaks != 0))
    else:
        ind=np.arange(len(c))
        
    #quiver plot
    
    py.clf()
    py.quiver(pts1[0,ind,0],pts1[0,ind,1],-diffx[ind],-diffy[ind])
    py.xlabel("x (mm)")
    py.ylabel("y (mm)")
    py.title("Distortion")
    #py.show()
    py.savefig("/Users/karr/"+prefix+"_distortion.jpg")

    #t1=np.partition(c.flatten(), -2)[-2]
    #t2=np.partition(c1.flatten(), -2)[-2]
    
    #plot map in mm
    py.clf()
    py.scatter(pts2[0,ind,0],pts2[0,ind,1],marker='s',c=c[ind],cmap='plasma')
    py.colorbar()
    py.title("Distortion (mm)")
    py.xlabel("x (mm)")
    py.ylabel("y (mm)")
    py.savefig("/Users/karr/"+prefix+"_distortion_col.jpg")

    #plot map in %
    py.clf()
    py.scatter(pts2[0,ind,0],pts2[0,ind,1],marker='s',c=c1[ind],cmap='plasma')
    py.colorbar()
    py.title("Distortion (% of field size)")
    py.xlabel("x (mm)")
    py.ylabel("y (mm)")
    py.savefig("/Users/karr/"+prefix+"_distortion_col1.jpg")

    py.clf()

    #t1=np.partition(fxs[ind].flatten(), -2)[-2]
    #t2=np.partition(fys[ind].flatten(), -2)[-2]
  

def plotVal(xs,ys,val,limit,plotrange,titl,prefix,suffix):

    if(limit > 0):
        #a quick cludge to filter out bad quality points in poorly focussed images
    
        ind=np.where((val < limit) & (val > 0) & (val > 0))
    else:
        ind=np.arange(len(c))

    #plot FWHMs

    py.clf()
    if(plotrange != None):
        sc=py.scatter(xs[ind],ys[ind],c=val[ind],marker="s",cmap='Purples',lw=0,s=20,vmin=plotrange[0],vmax=plotrange[1])

    else:
        sc=py.scatter(xs[ind],ys[ind],c=val[ind],marker="s",cmap='Purples',lw=0,s=20)

    py.colorbar(sc)
    py.title(titl)
    py.xlabel("X (mm)")
    py.ylabel("Y (mm)")
    py.savefig("/Users/karr/"+prefix+suffix+".jpg")

    
def to_fits(filename):

    #quick routine to convert raw to FITS

    a=np.fromfile(filename,dtype='uint16')
    image=a.reshape([5778,8960])

    pf.writeto(filename+".fits",image)
    print(image.min(),image.mean(),image.max())


def total_flux(xs,ys,image,x1,y1,scale):

    #back to pixel coordinates
    xc=xs*scale+x1
    yc=ys*scale+y1

    #simple photometry for total flux
    
    flux=total_flux(xc,yc,image,x1,y1,scale)

    
    ind=np.where(image > 0)
    
    #get average background
    back,a,b=sigmaclip(np.ravel(image[ind]))
    back=back.mean()

    
    flux=np.zeros((len(xc)))

    #cycle through the points, do a very simple aperture photometry
    #summing values in box
    
    for i in range(len(xc)):
    #for i in range(2):
        xx=np.round(xc[i])
        yy=np.round(yc[i])
        j1=int(xx)
        i1=int(yy)
        #print("circle point ",j1,i1)

        for ii in range(-4,5):
            for jj in range(-4,5):
                j1=int(xx+ii)
                i1=int(yy+jj)
                flux[i]=flux[i]+image[i1,j1]-back
        
    return flux


def mask_image(infile,outfile,x1,y1,x2,y2):

    image=pf.getdata(infile)
    image[0:x1,:]=0
    image[x2:,:]=0
    image[:,0:y1]=0
    image[:,y2:]=0
    pf.writeto(outfile,image,clobber=True)

