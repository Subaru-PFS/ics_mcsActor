from astropy.io.fits import getdata
import numpy as np

def getRegions(image,thresh1,thresh2,boxsize,xr,yr,xsize,ysize,nmin,nmax):

    print('Getting Regions\n')
    mask=np.zeros((xsize,ysize))
    xx=[]
    yy=[]
    rsize=[]
    x2=[]
    y2=[]
    xy=[]
    peak=[]
    for i in range(xr[0],xr[1]):
        for j in range(yr[0],yr[1]):

            #finds a point 
            if((image[i,j] > thresh1) and (mask[i,j]==0)):
                mask[i,j]=5
                #print(i,j)

                bx=0
                by=0
                bx2=0
                by2=0
                bxy=0
                npt=0
                tx=0
                ty=0
                #xb=0
                #yb=0
                maxval=thresh1;
                for ii in range(-boxsize,boxsize):
                    for jj in range(-boxsize,boxsize):
                        ip=i+ii
                        jp=j+jj
                        mask[ip,jp]+=1
                        #print('    ',ip,jp,i,j,ii,jj)

                        if(image[ip,jp] > thresh2):
                            #print('here',i,j,ip,jp)
                            mask[ip,jp]+=1
                            npt=npt+1
                            bx=bx+image[ip,jp]*ip
                            tx=tx+image[ip,jp]
                            by=by+image[ip,jp]*jp
                            ty=ty+image[ip,jp]
                            bx2=bx2+image[ip,jp]*ip*ip
                            by2=by2+image[ip,jp]*jp*jp
                            bxy=bxy+image[ip,jp]*ip*jp
                            if(image[ip,jp] > maxval):
                                maxval=image[ip,jp]
                            
                            #xb=xb+ip
                            #yb=yb+jp

                if(npt > nmin):
                    xx.append(bx/tx)
                    yy.append(by/ty)
                    x2.append(bx2/tx-(bx/tx)**2)
                    y2.append(by2/ty-(by/ty)**2)
                    xy.append(bxy/tx-(bx/tx)*(by/ty))
                    peak.append(maxval)
                    rsize.append(npt)
    print(str(len(xx))+' regions found\n')
    return xx,yy,x2,y2,xy,peak,rsize,mask
           
def callGetRegions(image):
    
    sz=image.shape
    xsize=sz[0]
    ysize=sz[1]
    thresh1=3000
    thresh2=1200
    nmin=10
    nmax=90
    
    xr=[120,4150]
    yr=[2660,6750]


    xr=[0,xsize]
    yr=[0,ysize]
    xx,yy,x2,y2,rsize,mask=getRegions(image,thresh1,thresh2,boxsize,xr,yr,xsize,ysize,nmin,nmax)

    return np.array(xx),np.array(yy),np.array(x2),np.array(y2),np.array(rsize),np.array(mask)

def windowedPos(image,x,y,boxsize,fwhmx,fwhmy,maxIt):
    
    threshold=2500

    #image - image array
    #xx,yy - positions returned by getRegion
    #swin d50/sqrt(8(log(2))
    
    diff=1e-4

    swinx=fwhmx/np.sqrt(8*np.log(2))
    swiny=fwhmy/np.sqrt(8*np.log(2))
    
    xwin=x
    ywin=y

    xmin=np.int(np.round(x-boxsize))
    ymin=np.int(np.round(y-boxsize))
    xmax=np.int(np.round(x+boxsize+1))
    ymax=np.int(np.round(y+boxsize+1))

    nIt=0
    dx=10
    dy=10
    while((dx > diff) & (dy > diff) & (nIt < maxIt)):

        xsum=0
        ysum=0
        nsumx=0
        nsumy=0
        
        #iterate over radius
        for i in range(xmin,xmax):
            for j in range(ymin,ymax):
                ri=np.sqrt((i-xwin)*(i-xwin)+(j-ywin)*(j-ywin))
                if(ri <= boxsize):
                    #if(image[i,j] > threshold):
                    #    print(ri,image[i,j])
                    wix=np.exp(-ri*ri/(2*swinx*swinx))
                    wiy=np.exp(-ri*ri/(2*swiny*swiny))
                    xsum+=wix*image[i,j]*(i-xwin)
                    ysum+=wiy*image[i,j]*(j-ywin)
                    nsumx+=wix*image[i,j]
                    nsumy+=wiy*image[i,j]
        xwin1=xwin+2*xsum/nsumx
        ywin1=ywin+2*ysum/nsumy
        #print(xsum,ysum,nsum,xwin1,ywin1)
        #print(xwin1-xwin,ywin1-ywin)

        dx=abs(xwin-xwin1)
        dy=abs(ywin-ywin1)
        nIt=nIt+1
        
        xwin=xwin1
        ywin=ywin1

    return xwin,ywin,nIt

def testCentroid(fwhmx,fwhmy,boxsize,xx,yy,rsize,mask,image,cfile):

    xval=[]
    yval=[]
    nint=[]
    #boxsize=4
    #fwhm=4
    #iterations=10
    maxIt=20

    for x,y in zip(xx,yy):
        
        xwin,ywin,nn=windowedPos(image,x,y,boxsize,fwhmx,fwhmy,maxIt)
        xval.append(xwin)
        yval.append(ywin)
        nint.append(nn)
    
    xval=np.array(xval)
    yval=np.array(yval)
    nint=np.array(nint)
    
    #centroidFile="../SEX/see_6077_6126_centroids.dat"
    centroidFile=cfile
    ff=np.loadtxt(centroidFile)

    ind=np.where(ff[:,0]==0)
    x1=ff[ind,2]
    y1=ff[ind,1]

    dx=[]
    dy=[]
    for x2,y2 in zip(xval,yval):
        dd=np.sqrt((x1-x2)**2+(y1-y2)**2)
        ind=np.where(dd==dd.min())
        dx.append((x2-x1[ind])[0])
        dy.append((y2-y1[ind])[0])

    return np.array(dx),np.array(dy),xval,yval,nint



def doCentroid(fname,boxsize,fwhmx,fwhmy,thresh1,thresh2,nmin,nmax,ff,filenum):
    
    boxsize1=7

    #read image
    print('reading image')
    image=getdata(fname)
    
    sz=image.shape
    xsize=sz[0]
    ysize=sz[1]
    xr=[0,xsize]
    yr=[0,ysize]

    #get regions
    xx,yy,x2,y2,xy,peak,rsize,mask=getRegions(image,thresh1,thresh2,boxsize1,xr,yr,xsize,ysize,nmin,nmax)

    xval=[]
    yval=[]
    nint=[]
    maxIt=20
    
    #do the centroiding and write to a file
    ii=0
    for x,y in zip(xx,yy):
        xwin,ywin,nn=windowedPos(image,x,y,boxsize,fwhmx,fwhmy,maxIt)
        print(filenum,xwin,ywin,x2[ii],y2[ii],peak[ii],xy[ii],nn,file=ff)

        xval.append(xwin)
        yval.append(ywin)
        ii=ii+1

    return np.array(xval),np.array(yval),mask
