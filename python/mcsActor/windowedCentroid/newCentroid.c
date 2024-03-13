#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h> 
#include <fitsio.h>
#include "centroid.h"
#include "centroid_types.h"
 #include <stdbool.h>
 
// create new node
struct QNode* newNode(long k)
{
    struct QNode* temp
        = (struct QNode*)malloc(sizeof(struct QNode));
    temp->key = k;
    temp->next = NULL;
    return temp;
}
 
// create empty queue
struct Queue* createQueue()
{
    struct Queue* q
        = (struct Queue*)malloc(sizeof(struct Queue));
    q->front = q->rear = NULL;
    return q;
}
 
// add a eleemnt to queue
void enQueue(struct Queue* q, long k)
{
    // new node
    struct QNode* temp = newNode(k);
 
    // If queue is empty, then new node is front and rear
    // both
    if (q->rear == NULL) {
        q->front = q->rear = temp;
        return;
    }
 
    // Add the new node at the end of queue and change rear
    q->rear->next = temp;
    q->rear = temp;
}
 
// remove element from queue
void deQueue(struct Queue* q)
{
    // If queue is empty, return NULL.
    if (q->front == NULL)
        return;
 
    // Store previous front and move front one node ahead
    struct QNode* temp = q->front;
 
    q->front = q->front->next;
 
    // If front becomes NULL, then change rear also as NULL
    if (q->front == NULL)
        q->rear = NULL;
 
    free(temp);
}


bool isValid(int *image, int mValue, long xsize, long ysize, long i, long j, long thresh)
{

  // is the pixel inside the region, above the threshold and unchecked in the mask

  //printf("%ld %ld %ld %ld %d %d %d\n",i,j,xsize, ysize, image[getInd2D(i,j,xsize)],mValue,thresh);
  if(i < 0 || i > xsize || image[getInd2D(i,j,xsize)] < thresh || mValue > 0)
    {
      return false;
    }
  else
    {

      return true;
    }
  
}

int running(int *image, long pixInd, long i, long j, long *bx, long *by, long *tt, long *bx2, long *by2, long *bxy, long *peak)
{

  (*bx) = (*bx) + image[pixInd] * i;
  (*by) = (*by) + image[pixInd] * j;
  (*tt) = (*tt) + image[pixInd];
  (*bx2) = (*bx2) + image[pixInd] * i * i;
  (*by2) = (*by2) + image[pixInd] * j * j;
  (*bxy) = (*bxy) + image[pixInd] * i * j;
  if(image[pixInd] > (*peak))
    {
      (*peak) = image[pixInd];
    }
  return 0;
}


struct cand_point *idSpot(int *image, int **mask, long xsize, long ysize, long i, long j, int thresh, int nmin)
{

  long bx=0;
  long by=0;
  long tt=0;
  long bx2=0;
  long by2=0;
  long bxy=0;
  long peak=0;

  long pixInd;
  long posX, posY;
  double npt = 0;
  struct Queue* xList = createQueue();
  struct Queue* yList = createQueue();
  struct cand_point *cand_curr;
  char *filename = "/Users/karr/dump4.txt";
  FILE *fp = fopen(filename, "w");

  enQueue(xList, i);
  enQueue(yList, j);
  long xRef = i;
  long yRef = j;
  pixInd = getInd2D(i,j,xsize);
  
  (*mask)[pixInd] = 1;
  int val;

			 
  int listSize = 1;
  while(listSize > 0)
    {

      // get the next positions and remove from queue
      posX = (xList->front)->key;
      posY = (yList->front)->key;
      deQueue(xList);
      deQueue(yList);

      //iterate until the list is empty
      listSize -=1;

      // check the neighbouring pixels for valid pixels, add to queue and update mask if valid.
      // we check the up, down, left, right pixels
      
      pixInd = getInd2D(posX+1, posY,xsize);
      if(isValid(image,(*mask)[pixInd],xsize,ysize,posX + 1,posY,thresh))
	{
	  // update the mask
	  (*mask)[getInd2D(posX + 1,posY,xsize)] = 1;

	  // add to the queue and update list size
	  enQueue(xList, posX + 1);
	  enQueue(yList, posY);
	  listSize +=1;

	  // now update the running totals for the first and second moments and peaks. subtract the reference pixels
	  // for sane math. 
	  val = running(image,pixInd,posX+1-xRef,posY-yRef,&bx, &by, &tt, &bx2, &by2, &bxy, &peak);
	  npt += 1.;
	}
      pixInd = getInd2D(posX-1,posY,xsize);
      if(isValid(image,(*mask)[pixInd],xsize,ysize,posX - 1,posY,thresh))
	{
	  (*mask)[getInd2D(posX - 1,posY,xsize)] = 1;
	  enQueue(xList, posX - 1);
	  enQueue(yList, posY);
	  listSize +=1;
	  val = running(image,pixInd,posX-1-xRef,posY-yRef,&bx, &by, &tt, &bx2, &by2, &bxy, &peak);
	  npt += 1.;
	}
      pixInd = getInd2D(posX,posY+1,xsize);
      if(isValid(image,(*mask)[pixInd],xsize,ysize,posX,posY+1,thresh))
	{
	  (*mask)[getInd2D(posX,posY+1,xsize)] = 1;
	  enQueue(xList, posX);
	  enQueue(yList, posY+1);
	  listSize +=1;
	  val = running(image,pixInd,posX-xRef,posY+1-yRef,&bx, &by, &tt, &bx2, &by2, &bxy, &peak);
	  npt += 1.;
	}
      pixInd = getInd2D(posX,posY-1,xsize);
      if(isValid(image,(*mask)[pixInd],xsize,ysize,posX,posY-1,thresh))
	{
	  (*mask)[getInd2D(posX,posY-1,xsize)] = 1;
	  enQueue(xList, posX);
	  enQueue(yList, posY-1);
	  listSize +=1;
	  val = running(image,pixInd,posX-xRef,posY-1-yRef,&bx, &by, &tt, &bx2, &by2, &bxy, &peak);
	  npt += 1.;
	}
	
    }


  // assign the values to the structure if above the minimum number of pixels. re-add on reference pixel to position. 
  double xval=((double)bx/(double)tt)+i;
  double yval=((double)by/(double)tt)+j;
  
  if((npt >= nmin))
    {	      
      cand_curr=(struct cand_point*)malloc(sizeof(struct cand_point));
      cand_curr->x=xval;
      cand_curr->y=yval;
      cand_curr->x2=(double)bx2/(double)tt-(double)bx/(double)tt*((double)bx/(double)tt);
      cand_curr->y2=(double)by2/(double)tt-(double)by/(double)tt*((double)by/(double)tt);
      cand_curr->xy=(double)bxy/(double)tt-(double)bx/(double)tt*((double)by/(double)tt);
      cand_curr->peak=(double)peak;
      cand_curr->npt=npt;
      
      return cand_curr;
    }

  // if there is no valid spot, return NULL
   else
    {
    return NULL;
    }
}
	

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
    
struct cand_point *getRegions(int *image,int thresh1,int thresh2,int boxsize,int boxsize1,int xsize,int ysize,int nmin,int **mask,int *npoints, int verbose)
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
     nmin, minimum  size of regions
     mask: integer mask file showing location of regions; for debugging purposes. 
     verbose: flag for diagnostic output

     returns: 

     list of structures containing the centroids and other parameters
     
     (also npoints: number of points found, mask: mask array for diagnostic purposes)

   */

  long i,j,ii,jj,ip,jp;
  //long bx,by,tt,bx2,by2,bxy,peak;
  int pixInd;
  double npt,xval,yval;
  struct cand_point *cand_head = NULL;
  struct cand_point *cand_curr = NULL;
  char *filename = "/Users/karr/dump1.txt";
  char *filename2 = "/Users/karr/dump3.txt";
  FILE *fp = fopen(filename, "w");
  FILE *fp1 = fopen(filename2, "w");
  (*npoints)=0;

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
	  if((image[pixInd] > thresh1) && ((*mask)[pixInd]==0))
	    {

	      cand_curr = idSpot(image, mask,xsize, ysize, i, j, thresh1, nmin);

	      if(cand_curr != NULL)
		{
		  cand_curr->next=cand_head;
		  cand_head=cand_curr;
		  (*npoints)+=1;
		}

	    }
	}
    }
  fclose(fp);
  fclose(fp1);
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
    int xmax=round(x+boxsize+1);
    int ymin=round(y-boxsize);
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

    //printf("EE %lf %lf %lf %lf\n",x,xwin,y,ywin);
    //finally, assign the results and pass back. 
    double *result=malloc(sizeof(double)*3);

    result[0]=xwin1;
    result[1]=ywin1;
    result[2]=nIt;

    return result;
  }
    

  
  
