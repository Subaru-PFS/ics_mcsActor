

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h> 
#include "mpfit.h"
#include "centroid_types.h"
#include "fitsio.h"
#include "centroid.h"

//Definitions for parallizing the code (NTHREAD=# of cores)

#define XSPLIT  1 //# of subregions in X direction
#define YSPLIT  1 //# of subregions in Y direction
#define NTHREAD 1  //# of cores

//toggle screen output for debugging/testing 

 
//--------------------------------------------------//

//array of structures for the loop
struct thread_data thread_data_array[NTHREAD];

//--------------------------------------------------//

/* Quick routine to print the results for a particular mpfit - label, number of iterations
   and the final value of the parameters. for debugging only */

void* subimage_Thread(void *arg)

{

  //do the find algorithm on a subimage

  long i,j,ii,jj;// counters

  /*sharp and round variables to reject non points*/

  double sharplim[2];  //LIMITS FOR SHARP PARAMETER (FIND)
  double roundlim[2];  //LIMITS FOR ROUND PARAMETER (FIND)

  /*set the variables for the routine from the structure (threaded part)*/

  int n_x=((struct thread_data*)arg)->n_x;         //x dimension of image
  int n_y=((struct thread_data*)arg)->n_y;         //x dimension of image
  int *image=malloc(n_x*n_y*sizeof(int));          //image

  int *h;        //convolved image
  int *imagemask;

  h = malloc(n_x*n_y*sizeof(int));
  imagemask = calloc(n_x*n_y,sizeof(int));

  int npoint=0;  //number of points found (debugging only)
  int ismax;     //is the point a max?

 
  //list to contain the guesses

  struct guess_point *guess_head = NULL;
  struct guess_point *guess_curr;

//struct cand_point *cand_head = NULL;
  struct cand_point *cand_curr;
 
  image=((struct thread_data*)arg)->image;            

  int fpix0=((struct thread_data*)arg)->fpix0;    //first and last pixels
  int fpix1=((struct thread_data*)arg)->fpix1;   

  int hmin=((struct thread_data*)arg)->hmin;          //data threshold
  double fwhm=((struct thread_data*)arg)->fwhm;          //fwhm
  int boxsize=((struct thread_data*)arg)->boxsize;      //box size for mpfit
  int nbox=((struct thread_data*)arg)->nbox;            //boxsize for find
  int pixels=((struct thread_data*)arg)->pixels;        //# of useable pixels in region (for find)
  int *mask=((struct thread_data*)arg)->mask;           //mask fo useable pixels in find
  double *gx=((struct thread_data*)arg)->gx;             //1d gaussian convolution kernel
  double *c1=((struct thread_data*)arg)->c1;             //weighted 2d kernel for round
  int nhalf=((struct thread_data*)arg)->nhalf;          //midpoint of image
  sharplim[0]=((struct thread_data*)arg)->sharplim[0];  //sharp and round limits
  sharplim[1]=((struct thread_data*)arg)->sharplim[1];
  roundlim[0]=((struct thread_data*)arg)->roundlim[0];
  roundlim[1]=((struct thread_data*)arg)->roundlim[1];
  struct cand_point *cand_head=((struct thread_data*)arg)->cand_list;  //output list of found points
  int fittype=((struct thread_data*)arg)->fittype;
  int verbose=((struct thread_data*)arg)->verbose;

  double sharp, around;
  int ix,iy;
  int hpoint=(boxsize-1)/2;


  if(verbose == 1)
    {
      printf("In Thread\n"); 
      printf("Starting Processing\n"); 
    printf("Convolving\n"); 
    }
  /*---------------------------------------------------------------*/

  
  //set up the image mask

  //make the mask and do the convolution
  imagemask=make_mask(image,nbox,n_x,n_y,hmin,verbose);
  h=convol_sep(image,imagemask,gx,n_x,n_y,nbox,verbose);

  if(verbose == 1)
  
    {
    printf("Finding Maxima\n"); 
    printf("hmin=%ld\n\n",hmin); 
    }

  /*---------------------------------------------------------------*/
  
  /*now for the find portion of the routine*/




  //now find pixels where value is greater than hmin and hmin is a local maximum
  for (i = 0;i< n_y;i++)
    {
  
      for (j=0;j<n_x;j++)
  	{

  	  //is the pixel greater than the threshold
  	  if(h[i*n_x+j]>hmin)
  	    {

  	      ismax=1;
  	      //is it a local maximum - compare with surrounding pixels
  	      for (ii=-1;ii<=1;ii++)
  		{
  		  for (jj=-1;jj<=1;jj++)
  		    {

  		      if(h[(i+ii)*n_x+(j+jj)] > h[i*n_x+j])
  			{
  			  ismax=0;
  			}
  		    }
  		}
  	      //if it's a local maximum, add to the coordinate list
  	      if(ismax==1)
  		{

		  //add to the list
		  guess_curr=(struct guess_point*)malloc(sizeof(struct guess_point));
		  guess_curr->x=j;
		  guess_curr->y=i;
		  guess_curr->next=guess_head;
		  guess_head=guess_curr;

  		  npoint=npoint+1;
  		}
	      
  	    }
  	}

    }
  
  if(verbose == 1)
    {
      printf("Found %i Maxima\n\n",npoint); 
      printf("Starting Statistics \n"); 
    }

  //free memory

  free(imagemask);
  free(h);

  /*---------------------------------------------------------------*/

  /*Now cycle through the guesses. first do round/sharp checks a la daofind,
    if the point passes, then do the centroiding and add it to the list*/

  

  /*cycle through the possible detections. idl/fortran codes in
  original daofind uses lots of goto replaced with if statements*/
  
  guess_curr=guess_head;

  while(guess_curr!=NULL)
  {
    //get the center points of the possible detection, and its value
    ix=guess_curr->x;
    iy=guess_curr->y;
    guess_curr=guess_curr->next;
  
    //check for points near the edge of the region. 
    if(ix > boxsize/2. && ix < n_x-boxsize/2. && iy > boxsize/2. && iy < n_y-boxsize/2.)
  	{

	  sharp=calculate_sharp(image,mask,nhalf,ix,iy,n_x,n_y,pixels,nbox);
          //if it fits the sharpness criteria continue to roundness
	  if(sharp < sharplim[1] && sharp > sharplim[0])
	    {
	      
	      around=calculate_round(image,mask,nhalf,ix,iy,n_x,n_y,pixels,nbox,c1);
	      if(around > roundlim[1] || around < roundlim[0])
		{
		  //printf("%d %d %lf\n",ix,iy,around);
		}

	  //if the roundness passes, add to the list
	  if(around > roundlim[0] && around < roundlim[1])
	    {

	      //we've got a good candidate, now centroiding. 
	      //pick the appropriate fit methods

	      if(fittype==0)
		{
		  cand_curr=getfit_simple(image,hpoint,n_x,n_y,ix,iy,fwhm,boxsize);
		}
	      else if(fittype==1)
		{
		  cand_curr=getfit_1d(image,hpoint,n_x,n_y,ix,iy,fwhm,boxsize);
		}
	      else if (fittype==2)
		{
		  cand_curr=getfit_2d(image,hpoint,n_x,n_y,ix,iy,fwhm,boxsize);
		}
	      cand_curr->next=cand_head;
	      cand_head=cand_curr;

	  //-----------------------------------------------------------//
	  
	    }
	    }
	} 
  }
  if(verbose == 1)
    {
      printf("Finished Filtering/Centroiding\n");
      printf("Tidying Results\n");
    }

  //now assign the list back to the thread variable to pass it back
  ((struct thread_data*)arg)->cand_list=cand_head;

  //free memory
  free(image);

  //exit the thread properly
  pthread_exit(0);
}



  /*-------------------------------------------------------------------------------------*/

struct centroids *centroid(int *image, int n_x, int n_y, int hmin, double fwhm,int boxsize,int *np, int fittype, double sharpLow, double sharpHigh, double roundLow, double roundHigh, int verbose)
{

  /*main routine. parses input values, calls setup routine, divides up the threads, 
    runs the threads, combines and filters the outputs */

  struct centroids *output; 
   
  int nbox=0;       //boxsize for find

  int ii,jj,i,j;          //counter variabls
  //int fname_check=0,hmin_check=0,fwhm_check=0;  //flags to check input

  double sharplim[2];  //sharp limits for find
  double roundlim[2];  //round limits for fine

  //default values for sharp/roundlim
  sharplim[0]=sharpLow;
  sharplim[1]=sharpHigh;
  roundlim[0]=roundLow;
  roundlim[1]=roundHigh;
  

  double *gx;       // 1D filter kernel
  int pixels=0;     // number of useable pixels in kernel (for sharp/round)
  int *mask;        //mask for useable pixels (for sharp/ round)
  int nhalf;        //dimension of kernel
  double *c1;       //kernel integrated along one dimension

  nbox=2*floor(0.637*fwhm) + 1;   //boxsize for smoothing

  if(verbose == 1)
    {
      printf("Starting Program\n");
      printf("  nbox = %d\n",nbox);
      printf("  fwhm = %lf\n",fwhm);
      printf("  thresh = %d\n",hmin);
      printf("  boxsize = %d\n",boxsize);

    }

  /*------------------------------------------------------------*/


  //Arrays for threading

  int ret[NTHREAD];       //RETURN VALUES
  pthread_t pth[NTHREAD]; //PTRHEAD VARIABLE

  struct cand_point *cand_list[NTHREAD];

  //first and last pixel values
  int fpix0;
  int lpix0;
  int fpix1;
  int lpix1;
  int nx,ny;

  //Some various pointers to keep track of the list
  struct cand_point *curr_val=NULL;    //first value
  struct cand_point *curr_pre=NULL;    //pointer previous to curr_val
  struct cand_point *check_val=NULL;   //second value (to compare to first value)
  struct cand_point *check_pre=NULL;   //pointer previous to check_val
  struct cand_point *top_val=NULL;     //top of list

  struct cand_point *cand_val;  //List f points

  double x1,y1;  //positions of first value
  double x2,y2;  //positions of first value
  double p1,p2;   //fluxes of first and second avalue
  double rs;     //square o fdistance between x1,y1 and x2,y2

  double rmin=(fwhm*3);  //radius to look for duplicate/false points
  int deadfirst=0;    //flag that first node has been deleted, to keep track fo pointers
  int firstval=1;     //flag that we're looking at the first node, to keep track of pointers
  char filename[sizeof "file100.fits"];

  //double *xpos,  *ypos,  *peak, *back, *fx, *fy, *qual;
  //xpos=malloc(ii*sizeof(double));
  //ypos=malloc(ii*sizeof(double));
  //fx=malloc(ii*sizeof(double));
  //fy=malloc(ii*sizeof(double));
  //qual=malloc(ii*sizeof(double));
  //peak=malloc(ii*sizeof(double));
  //back=malloc(ii*sizeof(double));

  //Allocate memory for kernel information
  gx = malloc(nbox*sizeof(double));  
  c1 = malloc(nbox*sizeof(double));  
  mask = malloc(nbox*nbox*sizeof(int));

  //Call routine to initialize kernel and information

  input_values(&nbox,&gx,&c1,&mask,&pixels,&nhalf,fwhm,verbose);
  //Verbose debugging information
  if(verbose ==1 )
    {
      printf("A %d %d %lf\n",pixels,nbox,fwhm);
      for (i=0;i<nbox;i++)
	{
	  printf("B %lf ",gx[i]);
	}
      printf("\n");
      
      for (i=0;i<nbox;i++)
	{
	  printf("C %lf ",c1[i]);
	}
      printf("\n");
      
      
      for (i=0;i<nbox;i++)
	{
	  for (j=0;j<nbox;j++)
	    
	    {
	      printf("%d ",mask[i+nbox*j]);
	    }
  	  printf("\n");
	  
	}
    }

  /*------------------------------------------------------------*/


  //Cycle through the subimages and set up the threading variables
  int ind, ind1, ind2;
  //Split the image in X and Y directions
  for (ii=0;ii<XSPLIT;ii++)
    {
      for (jj=0;jj<YSPLIT;jj++)
	{

	  ind=ii+XSPLIT*jj;
	  
	  cand_list[ind]=NULL;      //Initialize list of candidate points

	  /*Calculate boundaries of subimages. Different cases for subrejions on the edge of the image,
	    or interior, to properly calculate the overlaps. */

	  //Edge of the whole image at -X side 
	  if(ii==0)
	    {
	      fpix0=1;
	    }
	  //Split between interior split on -X side
	  else
	    {
	      fpix0=ii*n_x/XSPLIT+1-boxsize*5;
	    }
	  //Edge of the whole image at +X side 
	  if(ii==XSPLIT-1)
	    {
	      lpix0=n_x-1;
	    }
	  //Split between interior split on +X side
	  else
	    {
	      lpix0=(ii+1)*n_x/XSPLIT+boxsize*5;
	    }
	  //Edge of the whole image at -Y side 
	  if(jj==0)
	    {
	      fpix1=1;

	    }
	  //Split between interior split on -Y side
	  else
	    {
	      fpix1=jj*n_y/YSPLIT+1-boxsize*5;

	    }
	  //Edge of the whole image at +Y side 
	  if(jj==YSPLIT-1)
	    {
	      lpix1=n_y-1;
	    }
	  //Split between interior split on +Y side
	  else
	    {
	      lpix1=(jj+1)*n_y/YSPLIT+boxsize*5;
	    }

	  //Calculate size of subimage 
	  nx=lpix0-fpix0+1;
	  ny=lpix1-fpix1+1;

	  //Set the dimentions and image size for that thread
	  thread_data_array[ind].n_x=nx;
	  thread_data_array[ind].n_y=ny;
	  thread_data_array[ind].image=malloc(nx*ny*sizeof(int));
	  
	  //Assign the image data 
	  for(i=0;i<ny;i++)
	    {
	      for(j=0;j<nx;j++)
		{
		  ind1=i*nx+j;
		  //!! something is going wrong here, in the second index - fixed!
		  ind2=(i+fpix1-1)*n_x+(j+fpix0-1);
		  thread_data_array[ind].image[ind1]=image[ind2];
		}
	    }


	  
	  thread_data_array[ind].fpix0=fpix0;  
	  thread_data_array[ind].fpix1=fpix1;  
								      
	  //Set the values that need to be passed (same variable names as defined above)

	  thread_data_array[ind].hmin=hmin;        
	  thread_data_array[ind].fwhm=fwhm;        
	  thread_data_array[ind].boxsize=boxsize;
	  thread_data_array[ind].nbox=nbox;
	  thread_data_array[ind].nhalf=nhalf;
	  thread_data_array[ind].mask=mask;
	  thread_data_array[ind].pixels=pixels;
	  thread_data_array[ind].gx=gx;
	  thread_data_array[ind].c1=c1;
	  thread_data_array[ind].sharplim[0]=sharplim[0];
	  thread_data_array[ind].sharplim[1]=sharplim[1];
	  thread_data_array[ind].roundlim[0]=roundlim[0];
	  thread_data_array[ind].roundlim[1]=roundlim[1];
	  thread_data_array[ind].cand_list=cand_list[ind];
	  thread_data_array[ind].verbose=verbose;
	  thread_data_array[ind].fittype=fittype;
	  

	  //Set up the individual threads
	  ret[ind]=pthread_create(&pth[ind],NULL,subimage_Thread,(void *) &thread_data_array[ind]);

	}
    }

  /*Join the threads. Note that the threads do not share resources, so do not need to lock. Program will pause
    until all threads are done. */

    for (ii=0;ii<NTHREAD;ii++)
      {
	pthread_join(pth[ii],NULL);
      }


    /*End of thread part. Now join the outputs of each thread, and transform to coordinates of the whole 
      image. */

    //First go through the lists and link up each segment. 
    int iii;
    for(ii=0;ii<NTHREAD-1;ii++)
      {

	iii=0;

	//Set to the first value in the first list.
	if(ii==0)
	  {
	cand_val=thread_data_array[ii].cand_list;
	  }

	if(cand_val != NULL)  //check for empty list
	  {

	    while(cand_val->next != NULL)
	      {

		if(verbose==1)
		  {
		    printf("BB %d %lf %lf\n",ii,cand_val->x,cand_val->y);
		  }

		//Add the offset for that subimage
		if(iii != 0)
		  {
		    cand_val->x=cand_val->x+thread_data_array[ii].fpix0-1;
		    cand_val->y=cand_val->y+thread_data_array[ii].fpix1-1;
		  }       

		cand_val=cand_val->next;
		iii=iii+1;

	      }
	    //And the last point in the list
	    cand_val->x=cand_val->x+thread_data_array[ii].fpix0-1;
	    cand_val->y=cand_val->y+thread_data_array[ii].fpix1-1;

	    //and get the next list
	    cand_val->next=thread_data_array[ii+1].cand_list;
	  }
      }
    /*and the same for the last segment (or a single segment in the unthreaded case, when
      the above loop is not executed)*/

    cand_val=thread_data_array[NTHREAD-1].cand_list;

    if(cand_val != NULL)  //Are there items in the list?
      {
	while(cand_val != NULL)
	  {

	    if(verbose==1)
	      {
		//printf("ZZ %lf %lf\n",cand_val->x,cand_val->y);
	      }
	    cand_val->x=cand_val->x+thread_data_array[NTHREAD-1].fpix0-1;
	    cand_val->y=cand_val->y+thread_data_array[NTHREAD-1].fpix1-1;
	    
	    cand_val=cand_val->next;
	    
	  }
      }

    //Now filter out points with badly failed fits (fqual > 5)

    top_val=thread_data_array[0].cand_list;
    curr_val=top_val;
    curr_pre=top_val;

    while(curr_val != NULL)
      {

	if(curr_val->qual > 5)
	  {

	    //Delete First Node
	    if(curr_val==top_val)
	      {
		top_val=top_val->next;
		
		curr_val=curr_val->next;
		curr_pre=curr_pre->next;
	      }
	    //delete non first node
	    else
	      
	      {
		curr_pre->next=curr_val->next;
		curr_val=curr_val->next;
	      }

	  }
	else
	  {
	    curr_pre=curr_pre->next;
	    curr_val=curr_val->next;
	  }

      }


    //Now filter out duplicate points

    /*Reset to the beginning and set pointers. val and pre point to the same node at this point,
      as there is no previous node*/
    
    //top_val=thread_data_array[0].cand_list;
    curr_val=top_val;
    curr_pre=top_val;

    
   /*Filter out extaneous results. If two points are with x pixels of
   each other they are compared. If they are within 1 pixel of each
   other, it's a duplicate point, and the first one is
   retained. Otherwise, the point with the highest value is taken as
   the true point. */

  while(curr_val != NULL)
    {

      deadfirst=0;  //Reset Flag

      //get firest set of values
      x1=curr_val->x;
      y1=curr_val->y;
      p1=curr_val->peak;

      /*set the pointers to the next value. in the case of a single
	element list, both will be set to null*/
      check_val=curr_val->next;
      check_pre=curr_val;
      
      //now cycle through the comparisons
      while(check_val != NULL)
	{

	  //set the second set of values
	  x2=check_val->x;
	  y2=check_val->y;
	  p2=check_val->peak;

	  //Get the radius squared
	  rs=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
	  //if within the radius
	  if(rs < rmin)
	    {
	      /*We want to delete the second node - they are either the same point, duplicated,
		or the second point is fainter than the first*/

	      if((fabs(rs-rmin)<1) || (p1 > p2))  //delete pointed to node
		{

		  if(check_val->next != NULL)  //Not the last node
		    {
		      check_pre->next=check_val->next;
		      check_pre=check_val->next;
		      check_val=check_pre->next;
		    }
		  else  //case where it's the last node, so we don't try to follow a null pointer
		    {

		      check_pre->next=NULL;
		      check_val=check_pre->next;
		    }
 
		}
	      
	      //we want to delete the first node
	      else
		{
		  //special case when we're deleting the first element of the list
		  //set check_val to null to exit while loop
		  if(firstval==1)  
		    {

		      top_val=top_val->next;
		      check_val=NULL;
		      deadfirst=1;
		    }
		  else  
		    {

		      curr_pre->next=curr_val->next;
		      curr_val=curr_pre->next;

		      check_val=NULL;
		      deadfirst=1;
		    }
		}


	    }
	  //outside the radius, just implement points
	  else
	    {

	      check_val=check_val->next;
	      check_pre=check_pre->next;

	    }

	}

      
      //now we've compared an item, adn need to advance curr_val and curr_pre

      //we weren't on the first node
       if(firstval==0)  
	 {

	   //and we hadn't deleted the first node
	   if(deadfirst==0)
	     {	   

	       curr_pre=curr_pre->next;
	       curr_val=curr_val->next;
	     }
	 }
       /*special case for deleting the first node, to handle
	 increment of cand_val and cand_pre properly*/

       else
	 {

	   curr_val=curr_val->next;
	   //we didn't delete the first node, just curr_val needs to be incremented
	   if(deadfirst==0)  
	     {
	       firstval=0;
	     }
	   /*if we did, we need to increment both the pre and current
	     nodes, because the previous node is gone*/
	   else  
	     {
	       curr_pre=curr_pre->next;

	     }
	 }
    }


  //now it's sorted, print the results in region file format if we want to check it.

  if(verbose==1)
    {
      curr_val=top_val;
      while(curr_val != NULL)
	{
	  //printf("image;diamond point %lf %lf \n",curr_val->x,curr_val->y);
	  curr_val=curr_val->next;
	}

    }

  //count number of points
  ii=0;
  curr_val=top_val;
  while(curr_val != NULL)
    {
      curr_val=curr_val->next;
      ii=ii+1;
    }

  //assign the results to the appropriate variables to pass back
   

  output=malloc(sizeof(centroids)*ii);

  curr_val=top_val;
  for(i=0; i<ii;i++)
    {
      output[i].x=curr_val->x;
      output[i].y=curr_val->y;
      output[i].fx=curr_val->fx;
      output[i].fy=curr_val->fy;
      output[i].peak=curr_val->peak;
      output[i].back=curr_val->bg;
      output[i].qual=curr_val->qual;
     curr_val=curr_val->next;
    }

  np[0]=ii;


  return output;

}

/*--------------------------------------------------------*/
