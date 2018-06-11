
from matplotlib.backends.backend_pdf import PdfPages

import numpy as np
import pylab as py
import astropy.io.fits as pf


from centroid import centroid_only
import pickle

def check_values():

    times=['005','010','020','030','040','050','100']

    #get flat
    yvals=[1900,7100]
    xvals=[300,5400]

    flat=pf.getdata('flat_down_Master.fits')
    flat=flat[xvals[0]:xvals[1],yvals[0]:yvals[1]]
    flat=flat/flat.mean()

    fwhm=3.
    boxsize=9
    thresh=2500
    npoints=[]
    rrms=[]
    mmin=[]
    mmax=[]
    mmean=[]
    mmed=[]
    fname=[]

    n=1
    #basic stats
    for t in times:
        for i in range(n):
            fname="MCST_"+t+"_"+str(i+1).zfill(3)+".fits"

            data=pf.getdata(fname)
            
            py.clf()
            
            py.subplot(1,1,2)
            py.hist(data)
            data=data/flat
            py.subplot(2,1,2)

            py.imsave("hist_"+t+".jpg")

 
    #s=[]
    #
    #open("all_out.dat","w")
    #
    #pen(filename)
    #line in ff:
    #ff1=open("run2500_"+str(i+3).zfill(2)+".reg","w")
    #print(line)
    #data=pf.getdata(line.strip())
    #data=data[xvals[0]:xvals[1],yvals[0]:yvals[1]]
    #data=data/flat
    #
    #fname.append(line.strip())
    #mmin.append(data.min())
    #mmax.append(data.max())
    #mmean.append(data.mean())
    #mmed.append(np.median(data))
    #
    #rms=(data[np.where(data < data.max()/5)]).std()
    #rrms.append(rms)

        

def get_centroids():        
 
    cvals=[]

    ff2=open("all_out.dat","w")
    i=0
    ff=open(filename)
    for line in ff:
        ff1=open("run2500_"+str(i+3).zfill(2)+".reg","w")
        print(line)
        data=pf.getdata(line.strip())
        data=data[xvals[0]:xvals[1],yvals[0]:yvals[1]]
        data=data/flat
        
        fname.append(line.strip())
        mmin.append(data.min())
        mmax.append(data.max())
        mmean.append(data.mean())
        mmed.append(np.median(data))

        rms=(data[np.where(data < data.max()/5)]).std()
        rrms.append(rms)

        a=centroid_only(data.astype('<i4'),fwhm,thresh,boxsize)
        centroids=np.frombuffer(a,dtype='<f8')
        #reshape
        npoint=len(centroids)//7
        centroids=np.reshape(centroids,(npoint,7))

        print(centroids.shape)
        npoints.append(centroids.shape)

   
        cvals.append(centroids)
        for j in range(len(centroids)):
            print("circle point ",centroids[j][0]+yvals[0],centroids[j][1]+xvals[0],file=ff1)
            print(i,centroids[j][0]+yvals[0],centroids[j][1]+xvals[0],centroids[j][2],centroids[j][3],centroids[j][4],centroids[j][5],file=ff2)

        ff1.close()
        i=i+1
    ff2.close()

def get_rms():
    file1="all_out.dat"
    file2="all_rms.dat"

    ff=py.loadtxt(file1)

    
    #extract columns
    rnum=ff[:,0]
    x=ff[:,1]
    y=ff[:,2]
    fx=ff[:,3]
    fy=ff[:,4]
    peak=ff[:,5]
    back=ff[:,6]

    runs=int(rnum.max()+1)

    #extract the first exposure

    bind=np.where(rnum==0)
    xb=x[bind]
    yb=y[bind]
    fxb=fx[bind]
    fyb=fy[bind]
    peakb=peak[bind]
    backb=back[bind]

    nump=len(xb)
 
    #set up the arrays

    #distance from the first point

    dist=np.zeros((nump,runs))  

    #grids for the various values

    xall=np.zeros((nump,runs))
    yall=np.zeros((nump,runs))
    fxall=np.zeros((nump,runs))
    fyall=np.zeros((nump,runs))
    peakall=np.zeros((nump,runs))
    backall=np.zeros((nump,runs))

    #the rms for each point

    rms=np.zeros((nump))

    #set the first set of values to the first frame

    #now for each run
    for i in range(0,runs):

        #get the values for that run

        ind=np.where(rnum==i)
        xx=x[ind]
        yy=y[ind]
        fxx=fx[ind]
        fyy=fy[ind]
        peakk=peak[ind]
        backk=back[ind]
        
        nfound=0
        #now for each point in the reference frame
        for j in range(nump):
            found=0

            dd=np.sqrt((xb[j]-xx)*(xb[j]-xx)+(yb[j]-yy)*(yb[j]-yy))

            ind=np.where(dd==dd.min())
            
            if(dd.min() < 3):

                xall[j,i]=xx[ind]
                yall[j,i]=yy[ind]
                fxall[j,i]=fxx[ind]
                fyall[j,i]=fyy[ind]
                peakall[j,i]=peakk[ind]
                backall[j,i]=backk[ind]
    
    #dump the results to a file
    pickle.dump([xall,yall,fxall,fyall,peakall,backall],open(file2,"wb"))

def plot_rms():

    pref2="plots"
    pref3="1s A/C off"
    ##PLOTS - HISTORGRAM, LARGE MAP, SMALL MAPS (3), SAMPLE POS, SAMPLE FWHM
    ##MASTER PLOT - position of fiducial ma

    file2="all_rms.dat"

    #load data

    temp=pickle.load(open(file2,"rb"))
    xall=temp[0]
    yall=temp[1]
    fxall=temp[2]
    fyall=temp[3]
    peakall=temp[4]
    backall=temp[5]

    xb=xall[:,0]
    yb=yall[:,0]
    nruns=len(xall[0,:])
    runs=len(xall[0,:])

    #get rms
 
    nump=len(xall[:,0])

    rms=np.zeros((nump))
    nn=py.zeros((nump))
    fwhmx=np.zeros((nump))
    fwhmy=np.zeros((nump))

    for i in range(nump):
        ind=np.where((xall[i,:] > 0))
        dist=np.sqrt((xall[i,ind])**2+(yall[i,ind])**2)
        rms[i]=dist.std()

        fwhmx[i]=fxall[i,ind].mean()
        fwhmy[i]=fyall[i,ind].mean()
        nn[i]=dist.shape[1]
    
    #make plots


    with PdfPages(pref2+".pdf") as pdf:

        ##RMS colour coded by value
        py.figure()

        ind=np.where(rms > 0)
        py.hist(rms[ind],bins=0.01*np.arange(25))
        py.title(pref3+" - RMS")
        py.xlabel("X (pixels)")
        py.ylabel("RMS (pixels)")
        py.xlim([0,0.5])
        pdf.savefig()  
        py.close()
        
        py.figure()
        sc=py.scatter(xb,yb,c=rms,marker="s",cmap='Purples',lw=0,s=20,vmin=0,vmax=0.5)
        py.colorbar(sc)
        py.axis('equal')
        py.xlim([1900,7100])
        py.ylim([300,5400])
       
        py.title(pref3+" - RMS by Position")
        py.xlabel("X (pixels)")
        py.ylabel("Y (pixels)")
        py.savefig(pref2+".jpg")
        pdf.savefig()  
        py.close()

        py.figure()
    
        x1=np.array([1900,4490,3173])-1900
        x2=np.array([2650,5249,3983])-1900
        y1=np.array([3707,3687,1075])-300
        y2=np.array([4490,4497,1881])-300

        #for i in range(3):
        #    py.subplot(2,2,i+1)
        #    sc=py.scatter(xb,yb,c=rms,marker="s",cmap='Purples',lw=0,s=20,vmin=0,vmax=0.5)
        #    py.xlabel("X (pixels)",fontsize=8)
        #    py.ylabel("Y (pixels)",fontsize=8)
        #    py.title("Region "+str(i+1),fontsize=8)
        #    py.colorbar(sc)
        #    py.axis('equal')
        #
        #    py.xlim([x1[i]-1000,x2[i]-1000])
        #    py.ylim([y1[i]-300,y2[i]-300])
        #    py.tick_params(axis='x', labelsize=7)
        #    py.tick_params(axis='y', labelsize=7)
        # 
        #    
        #pdf.savefig()  


        py.figure()
        zz=fwhmx*fwhmy
        sc=py.scatter(xb,yb,c=zz,marker="s",cmap='Purples',lw=0,s=20,vmin=7,vmax=12)
        py.colorbar(sc)
        py.axis('equal')
        py.xlim([1900,7100])
        py.ylim([300,5400])

        #RMS by position

        py.title(pref3+" - FWHM(x*y) by Position")
        py.xlabel("X (pixels)")
        py.ylabel("Y (pixels)")
        pdf.savefig()  
        py.close()

        #py.figure()
        #ind1=np.where((nn < 29) & (nn > 25))
        #ind2=np.where(nn <= 25)

        ##dropouts
        #py.plot(xb[ind5],yb[ind5],'sb')
        ##py.plot(xb[ind2],yb[ind2],'dg')
        #py.title(pref3+" - Lower Matches")
        #py.xlabel("X (pixels)")
        #py.ylabel("Y (pixels)")
        #pdf.savefig()  
        #py.close()


        py.subplot(1,1,1)
        xx=np.zeros((nruns))
        yy=np.zeros((nruns))

        for i in range(nruns):
            x1=xall[:,i]
            y1=yall[:,i]
            
            ind=np.where(x1 > 0)

            xx[i]=(x1[ind]-xb[ind]).mean()
            yy[i]=(y1[ind]-yb[ind]).mean()

        py.plot(np.arange(nruns),xx,marker="d",markerfacecolor="g",linestyle="-",color="g",label="X")
        py.plot(np.arange(nruns),yy,marker="s",markerfacecolor="b",linestyle="-",color="b",label="Y")

        py.title(pref3+" - Average Position")
        py.xlabel("Frame Number")
        py.ylabel("P-P(0)")
        py.legend(framealpha=0.5)

        pdf.savefig()  
        py.close()

        #test images
        py.figure()

        xpos=[1000,1000,1000,2500,2500,2500,4500,4500,4500]
        ypos=[1000,2500,4500,1000,2500,4500,1000,2500,4500]
        
        #for i in range(len(xpos)):
        #    py.subplot(3,3,i+1)
        #    xp=xpos[i]
        #    yp=ypos[i]
        #    
        #    dd=np.sqrt((xb-xp)**2+(yb-yp)**2)
        #    ii=np.where(dd == dd.min())
        #    xx=xall[ii[0],:][0]
        #    yy=yall[ii[0],:][0]
        #    ind=np.where(xx > 0)
        #    py.plot(xx[ind]-xx[ind].mean(),'dg')
        #    py.plot(yy[ind]-yy[ind].mean(),'sb')
        #    py.tick_params(axis='x', labelsize=7)
        #    py.tick_params(axis='y', labelsize=7)
        #    py.title("("+str(xp)+","+str(yp)+")",fontsize=8)
        #py.suptitle(pref3+" - Position (Sample Fibres)")
        #py.annotate("Frame Number",xy=(0.4,0.03),rotation="horizontal",xycoords="figure fraction")
        #py.annotate("Position-mean(Position)",xy=(0.04,0.65),rotation="vertical",xycoords="figure fraction")
        #
        #pdf.savefig()  
        #py.close()
        

        py.subplot(1,1,1)
        ffx=np.zeros((nruns))
        ffy=np.zeros((nruns))

        for i in range(nruns):
            x1=fxall[:,i]
            y1=fyall[:,i]
            
            ind=np.where(x1 > 0)

            ffx[i]=(x1[ind]).mean()
            ffy[i]=(y1[ind]).mean()

        py.plot(np.arange(nruns-1)+1,ffx[1:],marker="d",markerfacecolor="g",linestyle="-",color="g",label="FWHM (X)")
        py.plot(np.arange(nruns-1)+1,ffy[1:],marker="s",markerfacecolor="b",linestyle="-",color="b",label="FWHM (Y)")

        py.title(pref3+" - Average FWHM")
        py.xlabel("Frame Number")
        py.ylabel("FWHM")
        py.legend(framealpha=0.5)

        pdf.savefig()  
        py.close()


        #for i in range(len(xpos)):
        #    py.subplot(3,3,i+1)
        #    xp=xpos[i]
        #    yp=ypos[i]
        #    dd=np.sqrt((xb-xp)**2+(yb-yp)**2)
        #    ii=np.where(dd == dd.min())
        #    zz=xall[ii[0],:][0]
        #    ind=np.where(zz > 0)
        #
        #    xx=fxall[ii[0],:][0]
        #    yy=fyall[ii[0],:][0]
        #
        #    py.plot(xx[ind]-xx[ind].mean(),'dg',label="X")
        #    py.plot(yy[ind]-yy[ind].mean(),'sb',label="Y")
        #    py.tick_params(axis='x', labelsize=7)
        #    py.tick_params(axis='y', labelsize=7)
        #    py.title("("+str(xp)+","+str(yp)+")",fontsize=8)
        #
        #py.suptitle(pref3+u" - FWHM (Sample Fibres)")
        #py.annotate("Frame Number",xy=(0.4,0.03),rotation="horizontal",xycoords="figure fraction")
        #py.annotate("FWHM-mean(FWHM)",xy=(0.04,0.65),rotation="vertical",xycoords="figure fraction")
        #
        #pdf.savefig()  
        #py.close()



        py.subplot(1,1,1)
        ffx=np.zeros((nruns))
        ffy=np.zeros((nruns))

        for i in range(nruns):
            x1=peakall[:,i]
            
            ind=np.where(x1 > 0)

            ffx[i]=(x1[ind]).mean()

        py.plot(np.arange(nruns-1)+1,ffx[1:],marker="d",markerfacecolor="b",linestyle="-",color="b",label="FWHM (X)")

        py.title(pref3+" - Peak Value Average")
        py.xlabel("Frame Number")
        py.ylabel("Peak")
        py.legend(framealpha=0.5)

        pdf.savefig()  
        py.close()

        
        py.figure()
        py.plot(xb,yb,'.k')
        py.plot(xpos,ypos,'or')
        py.text(1000,4000,'Region 1',color="white")
        py.text(3600,4000,'Region 2',color="white")
        py.text(2200,1300,'Region 3',color="white")
        py.title("Legend")
        py.savefig('legend.jpg')
        py.close()


#check_values('filelist.txt')
#get_rms()
plot_rms()
