import numpy as np
import cv2
import numpy.ma as ma


def approximateCoords(x,y,region):

    """

    Get lower left position and rotation angle from a set of pin-hole mask
    centroids.  This will likely break if the rotation angle is too high. 
    
    input

    x,y: 1D numpy arrays of coordinates (pixel)

    region: 4 element array with of the region of interest [x1,x2,y1,y2]

    returns: x1,y1 of lower left corner, angle in radians

    """

    #find points in the ergion
    
    ind=np.where((x > region[0]) & (x < region[1]) & (y > region[2]) & (y < region[3]))

    #find the minimum and maximum distance from the origin
    
    dd=x*x+y*y

    ind1=np.where(dd == dd.max())
    ind2=np.where(dd == dd.min())

    #calculate the angle of rotation wrt Y axis
    
    angle=np.arctan2((y[ind1]-y[ind2]),(x[ind1]-x[ind2]))-np.pi/4.

    return x[ind2],y[ind2],angle

def maskinMM(smallBox):

    """

    calculates the position of the holes in the pinhole mask, in mm
    
    set smallbox=1 to include the close spaced points, 0 to exclude

    returns x,y of holes

    note that this returns values with the lower left corner = (0,0), 
    subtract off the centre of (168,168) if needed
  
    """
    

    xx=[]
    yy=[]

    #region of close points
    
    pp1=np.array([6,13,8,15])+0.5
    pp3=np.array([28,35,8,15])+0.5
    pp2=np.array([17,24,30,37])+0.5

    #the dimensions of teh large spaced points are 8mm distance, 43 points
    
    for i in range(43):
        for j in range(43):
            iin=1
            ii=i
            jj=j

            #remove points in the close-spaced region
            
            if((ii > pp1[0]) & (ii < pp1[1]) & (jj > pp1[2]) & (jj < pp1[3])):
               iin=0
            if((ii > pp2[0]) & (ii < pp2[1]) & (jj > pp2[2]) & (jj < pp2[3])):
               iin=0
            if((ii > pp3[0]) & (ii < pp3[1]) & (jj > pp3[2]) & (jj < pp3[3])):

               iin=0
            if((i==21) and (j==21)):
                iin=0
            if(iin==1):
               #convert from point number to mm
                
               xx.append(ii*8)
               yy.append(jj*8)

    #same thing for the close spaced points if desired; 2 mm separation, 25 points
    
    if(smallBox==1):
        for i in range(25):
            for j in range(25):
                x1=i*2+(pp1[0]+0.5)*8
                y1=j*2+(pp1[2]+0.5)*8
                x2=i*2+(pp2[0]+0.5)*8
                y2=j*2+(pp2[2]+0.5)*8
                x3=i*2+(pp3[0]+0.5)*8
                y3=j*2+(pp3[2]+0.5)*8

                xx.append(x1)
                xx.append(x2)
                xx.append(x3)

                yy.append(y1)
                yy.append(y2)
                yy.append(y3)

    return np.array(xx),np.array(yy)

def scaleCentroids(x,y,x1,y1,scale):

    """

    scale the centroids to mm at mask, and shift to the origin

    input

    x1,y1: lower left spot
    x,y: coordinates of centroids
    scale: scale factor

    returns: transformed x,y

    """
    
    xc=(x-x1)/scale
    yc=(y-y1)/scale

    return xc,yc
 
def scaleMask(xx,yy,angle,flip):

    """

    rotate the mask as needed

    input

    xx,yy: hole positions
    angle: rotation angle
    flip: set to 1 to rotate by 90 degrees

    returns: transformed x,y

    """
    
    #apply any rotation (for matching purposes only)

    #add 90 degrees to rotation
    
    if(flip==1):
        angle=angle+np.pi/2

    #shift to centre, rotate, shift back
    xx=xx-168
    yy=yy-168
        
    xnew=xx*np.cos(angle)-yy*np.sin(angle)
    ynew=xx*np.sin(angle)+yy*np.cos(angle)
    
    xx=xnew+168
    yy=ynew+168

    return xx,yy


def matchPoints(xc,yc,xx,yy,fx,fy,peak):    

    """

    match the image spots to the mask holes

    input

    xc,yc: centroids
    xx,yy: hole coordinates
    fx,fy: fwhms
    peak: peak values

    returns: xs,ys,fxs,fyx,peaks: measured variables matched to the hole positions

    """

    #create arrays
    
    xs=np.zeros((len(xx)))
    ys=np.zeros((len(yy)))
    fxs=np.zeros((len(xx)))
    fys=np.zeros((len(xx)))
    peaks=np.zeros((len(xx)))
    
    #match points between two sets (nearest neighbour)
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



def getApproximateTransform(xx,yy):

    """

    Routine to take two sets of pinhole mask spots, and figure out the 
    transformation between them. Requires reasonably clean set of positions. 

    input: 
    xx,yy: spot positions

    output: 
    xm,ym: transformed mask coordinates
    xd1,yd1: x and y translations
    angle: rotation
    scale: scaling factor

    """

    #get mask coordinates, centred at origin
    x,y=maskinMM(1)
    x=x-168
    y=y-168
    
    #true mean values for mask

    mask_xm=x.mean()
    mask_ym=y.mean()

    #fine the maximum distance from the mask centre to a point (ie corner)
    dd=np.sqrt((x-mask_xm)**2+(y-mask_ym)**2)
    mask_s=dd.max()

    #get mean value of measured points
    
    cent_xm=xx.mean()
    cent_ym=yy.mean()

    #find the point farthest from the mean for the measured points (ie, corner)
    
    dd=np.sqrt((xx-cent_xm)**2+(yy-cent_ym)**2)
    ind1=np.where(dd==dd.max())
    cent_s=dd.max()

    #there are two equal values in the exact case, pick one at random
    
    x1=xx[ind1[0][0]]
    y1=yy[ind1[0][0]]

    #estimate the translaiton by subtracting the two sets of mean values
    
    xd1=cent_xm-mask_xm
    yd1=cent_ym-mask_ym

    #use the corner points to estimate the rotation. This has a 4 quadrant
    #uncertainty because we don't know the orientation yet. 
    
    cent_angle=np.arctan2(x1-cent_xm,y1-cent_ym)
    mask_angle=np.arctan2(x[0]-mask_xm,y[0]-mask_ym)

    angle=np.pi-cent_angle-mask_angle

    #get the scaling factor by comparing the two maximum distances

    scale=cent_s/mask_s

    #now transform the data to match the mask. 
    xa,ya=transformMask(xx,yy,-xd1,-yd1,0,1)
    x1,y1=transformMask(xa,ya,0,0,np.pi-angle,1/scale)

    #subtract any remaining translation using the minimum values of transformed
    #points
    
    xv=-x1.min()+x.min()
    yv=-y1.min()+y.min()
    
    x1=x1-xv
    y1=y1-yv


    #now we can break the degeneracy in quadrant using the fact that the
    #mask is not symmetric, and the average position in the y direction is
    #offset by 6.4 mm from the mask centre. 

    if(y1.mean() < -3):
        #we're good
        pass
    elif (y1.mean() > 3):
        x1,y1=transformMask(x1,y1,0,0,np.pi,1)
        angle=angle+np.pi
    elif (x1.mean() > 3):
        x1,y1=transformMask(x1,y1,0,0,-np.pi/2,1)
        angle=angle-np.pi/2
    elif (x1.mean() < -3):
        x1,y1=transformMask(x1,y1,0,0,np.pi/2,1)
        angle=angle+np.pi/2

    #now we're going to get the reverse transformation (mask to points)
    #so we can work in pixel coordinates from now on.
    
    xd1=xd1+xv*scale
    yd1=yd1+yv*scale

    xm,ym=transformMask(x,y,xd1,yd1,angle-np.pi,scale)

    return xm,ym,xd1,yd1,angle-np.pi,scale


def mask_image(infile,outfile,x1,y1,x2,y2):


    image=pf.getdata(infile)
    image[0:x1,:]=0
    image[x2:,:]=0
    image[:,0:y1]=0
    image[:,y2:]=0
    pf.writeto(outfile,image,clobber=True)


def transformMask(x,y,xd,yd,theta,s):

    """
    Apply a rigid transformation to the mask (trans, scale, rot). Mostly bookkeeping
    stuff. 

    input:
    x,y: mask positions
    xd,yd: translation
    theta: rotation (radians)
    s: scale

    output: transformed x,y

    """
    

    
    #create transformation matrix
    matrix=np.zeros((2,3))
    matrix[0,0]=np.cos(theta)*s
    matrix[0,1]=-np.sin(theta)*s
    matrix[1,0]=np.sin(theta)*s
    matrix[1,1]=np.cos(theta)*s
    matrix[0,2]=xd
    matrix[1,2]=yd

    #bookkeeping for coordinate format
    pts=np.zeros((1,len(x),2))

    pts[0,:,0]=x
    pts[0,:,1]=y

    pts=np.float32(pts)

    #do the transform
    pts1=cv2.transform(pts,matrix)

    #more bookkeeping
    xx=pts1[0,:,0]
    yy=pts1[0,:,1]

    return xx,yy


def matchAllPoints(centroidFile,xx,yy):

    """

    takes a batch of centroids created by getAllCentroids, registers
    the positions for each frame to the mask, and returns a set of arrays for
    each variable, in which points that were not detected have been masked. 

    xx, yy are the mask hole positions, and should be transformed to match the
    image

    input: 
    centroidFile: file with centroid output
    xx,yy: transformed mask coordinates

    output: 
    xArray,yArray: array of registered positions
    fxArray,fyArray,backArray,peakArray: registered parameters for spots

    """

    #create arrays
    centroids=np.loadtxt(centroidFile)

    #size of output array
    npoints=len(xx)
    nfiles=np.int(centroids[:,0].max())+1

    
    xArray=np.zeros((npoints,nfiles))
    yArray=np.zeros((npoints,nfiles))
    fxArray=np.zeros((npoints,nfiles))
    fyArray=np.zeros((npoints,nfiles))
    peakArray=np.zeros((npoints,nfiles))
    backArray=np.zeros((npoints,nfiles))

    #load the centroid data
    
    print(str(nfiles)+" frames. Matching ",end="")
    #for each image
    for i in range(nfiles):

        print(str(i+1)+", ",end="")

        #get spot positions, etc from a particular image

        ind=np.where(centroids[:,0]==i)

        ll=centroids[ind,1].shape[1]
        x=centroids[ind,1].reshape(ll)
        y=centroids[ind,2].reshape(ll)
        fx=centroids[ind,3].reshape(ll)
        fy=centroids[ind,4].reshape(ll)
        peak=centroids[ind,5].reshape(ll)
        back=centroids[ind,6].reshape(ll)


        #nearest neighbour matching
        for j in range(npoints):

            dd=np.sqrt((xx[j]-x)*(xx[j]-x)+(yy[j]-y)*(yy[j]-y))
            ind=np.where(dd==dd.min())

            #need to filter in case points are missing
            if(dd.min() < 30):

                xArray[j,i]=x[ind]
                yArray[j,i]=y[ind]
                fxArray[j,i]=fx[ind]
                fyArray[j,i]=fy[ind]
                peakArray[j,i]=peak[ind]
                backArray[j,i]=back[ind]
    print()
    #mask unfound values
    xArray=ma.masked_where(peakArray == 0,xArray)
    yArray=ma.masked_where(peakArray == 0,yArray)
    fxArray=ma.masked_where(peakArray == 0,fxArray)
    fyArray=ma.masked_where(peakArray == 0,fyArray)
    backArray=ma.masked_where(peakArray == 0,backArray)
    peakArray=ma.masked_where(peakArray == 0,peakArray)
    
    return xArray,yArray,fxArray,fyArray,backArray,peakArray


def getTransform(x,y,xx,yy,getVals):

    """

    given two sets of registered points, estimate the rigid transformation
    this is in a separate routine mostly for bookkeeping purposes. 
    Returns transformation matrix, and if getVals == 1 returns the 
    extracted parameters (rotation, translation, scale) as well. 

    input:
    x,y: input positions
    xx,yy: transformed positions
    getVales: if ==1, return parameters too

    output: 
    transformation: matrix 
    xd,yd: translations
    sx,sy: scalings
    rotation: rotation (radians)

    """
    
    pts1=np.zeros((1,len(x),2))
    pts2=np.zeros((1,len(x),2))

    pts1[0,:,0]=x
    pts1[0,:,1]=y

    pts2[0,:,0]=xx
    pts2[0,:,1]=yy

    
    pts1=np.float32(pts1)
    pts2=np.float32(pts2)

    transformation = cv2.estimateRigidTransform(pts1, pts2, False)

    if(getVals == 0):
        return transformation
    if(getVals == 1):


        sx=np.sqrt(transformation[0,0]**2+transformation[0,1]**2)
        sy=np.sqrt(transformation[1,0]**2+transformation[1,1]**2)

        xd=transformation[0,2]
        yd=transformation[1,2]

        rotation = np.arctan2(transformation[1,0]/sx,transformation[1,1]/sy)
        
        return transformation,xd,yd,sx,sy,rotation

