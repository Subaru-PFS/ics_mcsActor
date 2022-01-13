#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h> 
#include "fitsio.h"
#include "centroid.h"
#include "centroid_types.h"


int *getParams(struct cand_point *cand_list,double *fwhmx,double *fwhmy)
{
  double n=0;

  struct cand_point *cand_curr=cand_list;
  while(cand_curr!=NULL)
    {
      (*fwhmx)+=cand_curr->x2;
      (*fwhmy)+=cand_curr->y2;
      n+=1;
      cand_curr=cand_curr->next;
    }
  (*fwhmx)=(*fwhmx)/n;
  (*fwhmy)=(*fwhmy)/n;
	  
  return 0;

}


int getInd2D(int j,int i,int size)
  {

    //index management for 2D arrays
    return i*size+j;
  }


int maxValI(int val1,int val2)
{

  /*return maximum of two values (integer) */
  
  if(val1 > val2)
    {
      return val1;
    }
  else
    {
      return val2;
    }
}


double maxValD(double val1,double val2)
{

  /*return maximum of two values (double) */

  if(val1 > val2)
    {
      return val1;
    }
  else
    {
      return val2;
    }
}
    
struct cand_point *getRegions(int *image,int thresh1,int thresh2,int boxsize,int boxsize1,int xsize,int ysize,int nmin,int *mask,int *npoints, int verbose)
{

  /* when passed an image, finds regions above the threshold. The
     assumption is that the regions are isolated in pixel space, and
     therefore blending is not an issue. As a result, the contiguity of the 
     pixels is not explicitly checked for; it is assumed that pixels in the PSF will 
     be within 2*boxsize pixels of the first detected pixel.

     input:

     image - input image;
     xsize, ysize  dimensions of the image
     thresh1 - threshold for finding a pixel to start the region
     thresh2 - threshold for finding pixels in the region 
             - these two are different in order to start 
	       the detection close to the centre of the PSF, rather than at the edge.

     boxsize: boxsize for searching for pixels. 
     nmin, minimum and maximum size of regions
     mask: integer mask file showing location of regions; for debugging purposes. 
     verbose: flag for diagnostic output

     returns: 

     list of structures containing the centroids and other parameters
     
     (also npoints: number of points found, mask: mask array for diagnostic purposes)

   */

  long i,j,ii,jj,ip,jp;
  long bx,by,tt,bx2,by2,xb,yb,bxy,peak;
  int pixInd;
  double npt,xval,yval;
  struct cand_point *cand_head = NULL;
  struct cand_point *cand_curr;

  (*npoints)=0;

  //printf("VV %d %d %d %d %d %d %d %d\n",thresh1,thresh2,xsize,ysize,boxsize,nmin,verbose);
  
  if(verbose==1)
    {
      printf("Finding Regions");
    }
  //cycle through the image

  for (i=0;i<xsize;i++)
    {
      for (j=0;j<ysize;j++)
	{
	  //first, find a pixel in the region, using the higher threshold
	  pixInd=getInd2D(i,j,xsize);
	  //is it above the threshold, and not previously chedked
	  if((image[pixInd] > thresh1) && (mask[pixInd]==0))
	    {

	      //initialize variables
	      npt=0.;
	      mask[pixInd]=5;
	      bx=0;
	      by=0;
	      tt=0;
	      bx2=0;
	      by2=0;
	      bxy=0;
	      peak=thresh1;
	      
	      //now cycle around the detcted point using the lower threshold.
	      for(ii=-boxsize;ii<=boxsize;ii++)
		{
		  for(jj=-boxsize;jj<=boxsize;jj++)
		    {
		      
		      ip=ii+i;
		      jp=jj+j;
		      //check for edges of image
		      if((ip >= 0) && (ip < xsize) && (jp >= 0) && (jp < ysize))
			{

			  //update mask
			  pixInd=getInd2D(ip,jp,xsize);
			  mask[pixInd]+=1;
			 
			  //now check for the threshold
			  if(image[pixInd] >= thresh2)
			    {

			      mask[pixInd]+=1;

			      //keep a running tally of for the isophotal
			      //centroids and shape parameters for each region
			      npt+=1.;
			      bx+=image[pixInd]*ii;
			      by+=image[pixInd]*jj;
			      tt+=image[pixInd];
			      bx2+=image[pixInd]*ii*ii;
			      by2+=image[pixInd]*jj*jj;
			      bxy+=image[pixInd]*ii*jj;

			      //calculate the peak value
			      if(image[pixInd] > peak)
				{
				  peak=image[pixInd];
				}

			    }
			}
		      
		    }
		}

	      xval=((double)bx/(double)tt)+i;
	      yval=((double)by/(double)tt)+j;
	      //now that we have a region, check its size to filter out hot pixels etc.
	      if((npt >= nmin) && ((xval-boxsize1) > 0) && ((yval-boxsize1) > 0) && ((xval + boxsize1) < xsize) &&  ((yval + boxsize1) < ysize))
		{

		  //assign to the structure
		  cand_curr=(struct cand_point*)malloc(sizeof(struct cand_point));

		  cand_curr->x=xval;
		  cand_curr->y=yval;
		  cand_curr->x2=((double)bx2/(double)tt-((double)bx/(double)tt)*((double)bx/(double)tt));
		  cand_curr->y2=((double)by2/(double)tt-((double)by/(double)tt)*((double)by/(double)tt));
		  cand_curr->xy=(double)bxy/(double)(tt)-((double)bx/(double)tt)*((double)by/(double)tt);
		  ////printf("%d %ld %ld %ld %ld %ld %ld %lf %lf %lf %lf \n",(*npoints),bx,by,tt,bx2,by2,bxy,cand_curr->x,cand_curr->y,cand_curr->x2,cand_curr->y2);
		  cand_curr->qual=0;
		  cand_curr->peak=peak;
		  cand_curr->npt=npt;
		  cand_curr->next=cand_head;
		  cand_head=cand_curr;

      		  (*npoints)+=1;

		}
	    }
		
	}
    }
  if (verbose==1)
    {
      //printf("Found %d Regions",(*npoints));
    }
  return cand_head;

}

double *windowedPos(int *image,double x, double y,int boxsize,double fwhmx, double fwhmy,int maxIt, int xsize, int ysize,int verbose)
  {

    /*
      calculate windows posiitions for a single PSF region. Based on the SEXtractor
      windowed parameters (ref). 

      input
        image: image
	xsize,ysize: size of image
	x,y: initial guesses for central position (isophotal values)
	boxsize: boxsize for windowing
	fwhm: estimate of the FWHM. 
	maxIt: maximum number of iterations. 
	precision: required precision for iteration

     */
    //required precision
    double precision=1e-4;
    
    //parameter for windowing
    double swinx=fwhmx/sqrt(8*log(2));
    double swiny=fwhmy/sqrt(8*log(2));

    //initialize the variables
    
    double xwin=x;
    double ywin=y;

    int nIt=0;
    double dx=10;
    double dy=10;

    int pixInd;

    double xsum, ysum, nsumx, nsumy, ri, wix, wiy;
    int boxsize2=boxsize*boxsize;

    int i,j,ii,jj;
    
    int xmin=round(x-boxsize);
    int xmax=round(x+boxsize);
    int ymin=round(y-boxsize+1);
    int ymax=round(y+boxsize+1);

    double xwin1,ywin1;

    
    //cycle through until precision is met *or* maxIt is reached
    while((dx > precision) && (dy > precision) && (nIt < maxIt))
      {
	xsum=0;
	ysum=0;
	nsumx=0;
	nsumy=0;


	//sum over the region
	for(i=xmin;i<=xmax;i++)
	  {
	    for(j=ymin;j<ymax;j++)
	      {
		//circular aperture
		ri=(i-xwin)*(i-xwin)+(j-ywin)*(j-ywin);
		if(ri < boxsize2)
		  {
		    //calculate the values
		    pixInd=getInd2D(i,j,xsize);;
		    
		    wix=exp(-ri/(2*swinx*swinx));
		    wiy=exp(-ri/(2*swiny*swiny));
		    
		    xsum+=wix*image[pixInd]*(i-xwin);
		    ysum+=wiy*image[pixInd]*(j-ywin);
		    nsumx+=wix*image[pixInd];
		    nsumy+=wiy*image[pixInd];
		  }
	      }
	  }
	//new value, check for precision

	xwin1=xwin+2*xsum/nsumx;
	ywin1=ywin+2*ysum/nsumy;

	dx=fabs(xwin-xwin1);
	dy=fabs(ywin-ywin1);
	nIt=nIt+1;

	xwin=xwin1;
	ywin=ywin1;
      }

    //finally, assign the results and pass back. 
    double *result=malloc(sizeof(double)*3);

    result[0]=xwin1;
    result[1]=ywin1;
    result[2]=nIt;

    return result;
  }
    

  
  
