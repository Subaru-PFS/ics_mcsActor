"""

Centroiding related routines

"""


import numpy as np
import numpy.ma as ma
from scipy.stats import sigmaclip

import vis_util as visutil

import centroid as centroid

def getCentroids(image,thresh1,thresh2,fwhmx,fwhmy,boxFind,boxCent,nmin,nmax,maxIt):

    """
    
    wrapper to run centroid command, and convert the output to arrays

    input: 
    
    image (numpy array)
    thresh1, thresh2: thresholds for spot detection
    fwhmx,fwhmy: default 0, set to value ot override automatic determinatino
    boxFind, boxCent - box size for finding and centroids
    nmin, nmax - min and max accpetable # of pixels in spot
    maxIt: maximum number of iterations in centroiding


    output:

    x,y: positions
    fx, fy: FWHMs
    back: background
    peak: peak values


    """

    #centroid and reshape the results
    
    a=centroid.centroid_only(image.astype('<i4'),fwhmx,fwhmy,thresh1,thresh2,boxFind,boxCent,nmin,nmax,maxIt,0)

    centroids=np.frombuffer(a,dtype='<f8')
    centroids=np.reshape(centroids,(len(centroids)//7,7))

    #extract useful values

    x=centroids[:,0]
    y=centroids[:,1]
    fx=(centroids[:,2])
    fy=(centroids[:,3])
    back=(centroids[:,4])
    peak=(centroids[:,5])
    qual=(centroids[:,6])

    return x,y,fx,fy,back,peak,qual


def getAllCentroidsDB(conn,frameIDs):

    """

    retrieves a set of centroids from the database, for a sequence of frameIDs

    Input:  
       conn: database connection
       frameIDs: list of frame id numbers

    """
    
    #for commissioning with MCS
    
    moveId = 1
    
    #make a blank array for the centroid array
    centroids=np.array([])
    
    i=0
    
    #cycle through each ID number
    for id in frameIDs:
        
        #SQL for getting a set of centroids
        #cmd_string = f"""select * from mcsEngTable where frameId={id} and moveId=1"""???
        cmd_string=""
        data=np.array([]) 
        n = 0
        with conn.cursor() as curs:
                curs.execute(cmd_string)
                rows=curs.fetchall()
                for idx, val in enumerate(rows):
                    if idx == 0: data = val 
                    if idx != 0: data = np.vstack([data,val])
        conn.commit()
        
        
        #some data massaging into the right form. 
        cen=data[:,5:11]
        cen1=np.zeros((cen.shape[0],7))
        
        #add an index to the first number
        cen1[:,0]=i
        
        #copy over the centroids
        cen1[:,1:7]=cen
        
        #create master array
        if(i==0):
            centroids=cen1
        else:
            centroids=np.concatenate((centroids,cen1),axis=0)

    return centroids    


def getAllCentroids(files,outfile,thresh1,thresh2,fwhmx,fwhmy,boxFind,boxCent,nmin,nmax,maxIt):        

    """

    wrapper to do the centroiding for a set of files. Writes to a text
    output file.  Returns the last set of xy positions (used for the
    next routine.

    input: 
    
    files: list of FITS files
    outfile: name of output text file
    FWHM, boxsize, thresh, rl,rh,sl,sh: parameters for centroiding routine

    output: writes to a file with each line of the format
       framenum, x, y, fx, fy, peak, back, qual

    """

    #text file for output
    ff=open(outfile,"w")
    nfiles=len(files)
    
    #cycle through files
    filenum=0

    print(str(nfiles)+" Frames. Centroiding ",end="")
    for file in files:
        #read in image
        image=visutil.getImage(file)
        #centroid
        x,y,fx,fy,back,peak,qual=getCentroids(image,thresh1,thresh2,fwhmx,fwhmy,boxFind,boxCent,nmin,nmax,maxIt)
        print(str(filenum+1)+": "+str(len(x))+", ",end="")
    
        #write results to file, labelling the file number for bookkeeping
        for i in range(len(x)):
            print(filenum,x[i],y[i],fx[i],fy[i],back[i],peak[i],qual[i],file=ff)
        filenum+=1
            
    print()
    ff.close()
    return x,y


def getThreshT(image,sigma1,sigma2):

    #sigma clip the image to get rid of the spots
    subIm=image
    a,b,c=sigmaclip(subIm.ravel(),4,4)

    #sigma clipped mean + rms*factor
    thresh1=a.mean()+a.std()*70
    thresh2=a.mean()+a.std()*15

    return thresh1,thresh2

    
def getThresh(image,xrange,yrange,sigma1,sigma2):

    """

    calculate the two thresholds used by the centroiding routine.
    
    input:
       image: image
       xrange,yrange: upper and lower limit of the region of the mask
       sigma1: higher sigma multiplier for finding region
       sigma2: lower sigma multiplier for calculating weighted centroid

    returns:
       thresh1, thresh2: upper and lower thresholds for finding spots,
       in pixel flux values.
    
    """

    #sigma clip the image to get rid of the spots
    subIm=image[xrange[0]:xrange[1],yrange[0]:yrange[1]]
    a,b,c=sigmaclip(subIm.ravel(),4,4)

    #sigma clipped mean + rms*factor
    thresh1=a.mean()+a.std()*sigma1
    thresh2=a.mean()+a.std()*sigma2

    return thresh1,thresh2,xrange,yrange


def getRegion(image,high,factor):

    im,a,b=sigmaclip(image,high=high)
    xlow,xhigh,ylow,yhigh=getBoundary(image,a,b,0)
    im,a,b=sigmaclip(image[xlow:xhigh,ylow:yhigh],high=high)
    rrms=im.std()/factor
    xlow,xhigh,ylow,yhigh=getBoundary(image,a,b,rrms)

    return [xlow,xhigh],[ylow,yhigh]

def getBoundary(image,a,b,rms):

    """

    find the region that is occupied by the pinhole mask. 

    input:
 
        image
        a,b - upper and lower limits for sigma clipping
        rms - adjustment factor for short axis of image (default zero)

    output:

        x1,y1,x2,y2: limits of box for the pinhole mask structure. 

    """

    #set up the variables
    sz=image.shape

    prof1=np.zeros((sz[0]))
    prof0=np.zeros((sz[1]))

    
    #do a sigma clipped summing collapsing along each axis.
    #there's probabably a fancy python way to make this more
    #efficient
    
    for i in range(sz[0]):
        ind=np.where(image[i,:] < b)
        prof1[i]=np.mean(image[i,ind])
    for i in range(sz[1]):
        ind=np.where(image[:,i] < b)
        prof0[i]=np.mean(image[ind,i])

        
    #get the mean of each profile, with an adjustment for the short axis
    #of the image
    
    pm0=prof0.mean()
    pm1=prof1.mean()-rms

    #next step - move along the summed profile and find the step function
    #portion by looking for the place in which the step crossed the (adjusted)
    #mean value. 

    #start from the left and move right, then start from the right and move left.
    #if we reach the middle of the image without finding it, then we assume
    #the region is right up against the edge of the image. 
    
    ###########################################
    
    found=0
    i=0
    while((found==0) & (i < sz[0]/2)):
        if((prof1[i]-pm1)*(prof1[i+1]-pm1)<0):
            found=1
        i=i+1
        
    x1=i
        
    ###########################################3
    i=int(sz[0]-1)
    found=0
    while((found==0) & (i > sz[0]/2)):
        if((prof1[i]-pm1)*(prof1[i-1]-pm1)<0):
            found=1
        i=i-1

    if(found==1):
        x2=i
    else:
        x2=sz[0]-1
        
    ###########################################3

    #start at 1000 rather than zero in the long axis
    #to avoid variation in background flux at the edge of the image
    
    found=0
    i=1000
    
    while((found==0) & (i < sz[1]/2)):
        if((prof0[i]-pm0)*(prof0[i+1]-pm0)<0):
            found=1
        i=i+1
        
    y1=i
    
    ###########################################3
    i=sz[1]-1
    found=0
    while((found==0) & (i > 1)):
        if((prof0[i]-pm1)*(prof0[i-1]-pm0)<0):
            found=1
        i=i-1
        
    y2=i

    ###########################################3

    return x1,x2,y1,y2

def getAutoThresh(image,xr,yr,mx,my):

    #get centre of mask in mm at pfi (offset from centre)
    xyin=np.array([[mx-336/2,mx+336/2],[my-336/2,my+336/2]])

    #transform into pixels
    xyout=CoordTransp.CoordinateTransform(xyin,0,'pfi_mcs_wofe',inr=-angle*180/np.pi,cent=np.array([[xr],[yr]]))

    #centre of mask in pixels
    mpx=xyout[0,0]
    mpy=xyout[0,1]

    xx,yy=np.meshgrid(np.arange(image.shape[0]),np.arange(image.shape[1]))
    
    dfromc=np.sqrt((xx-mpx)**2+(yy-mpy)**2)

    ind=np.where(dfromc < 1500)

    a,b,c=sigmaclip(image[ind],4,4)

    thresh1=a.mean()+a.std()*70
    thresh2=a.mean()+a.std()*15

    return thresh1,thresh2


def dist(x1,x2,y1,y2):

    dd=np.sqrt((x1-x2)**2+(y1-y2)**2)
    return dd
