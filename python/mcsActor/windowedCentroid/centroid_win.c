

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h> 
#include "centroid_types.h"
#include "fitsio.h"
#include "centroid.h"

//Definitions for parallizing the code (NTHREAD=# of cores)

#define XSPLIT  3//# of subregions in X direction
#define YSPLIT  4//# of subregions in Y direction
#define NTHREAD 12//# of cores

//toggle screen output for debugging/testing 

 
//--------------------------------------------------//

//array of structures for the loop
struct thread_data thread_data_array[NTHREAD];


void freeAll(struct cand_point **head)
{

  /*routine to free the list of structures*/

  struct cand_point* current = *head;
  struct cand_point* next;
  while (current != NULL)
    {
      next = current->next;
      free(current);
      current=next;
    }
  *head=NULL;
}
//--------------------------------------------------//

/* Quick routine to print the results for a particular mpfit - label, number of iterations
   and the final value of the parameters. for debugging only */

void* subimage_Thread(void *arg)
{

  //do the find algorithm on a subimage

  long i,j,ii,jj;// counters
  /*set the variables for the routine from the structure (threaded part)*/

  int n_x=((struct thread_data*)arg)->n_x;         //x dimension of image
  int n_y=((struct thread_data*)arg)->n_y;         //x dimension of image
  int *image=malloc(n_x*n_y*sizeof(int));          //image

  int *imagemask;

  imagemask = calloc(n_x*n_y,sizeof(int));
  //printf("a  n_x = %d\n",n_x);
  //printf("a  n_y = %d\n",n_y);

  int npoint=0;  //number of points found (debugging only)
  int ismax;     //is the point a max?
 
  //list to contain the guesses

//struct cand_point *cand_head = NULL;
  struct cand_point *cand_curr;
 
  image=((struct thread_data*)arg)->image;            

  int fpix0=((struct thread_data*)arg)->fpix0;    //first and last pixels
  int fpix1=((struct thread_data*)arg)->fpix1;   

  int thresh1=((struct thread_data*)arg)->thresh1;          //data threshold
  double fwhmx=((struct thread_data*)arg)->fwhmx;          //fwhm
  int boxFind=((struct thread_data*)arg)->boxFind;      //box size for mpfit
  int nmin=((struct thread_data*)arg)->nmin;          //data threshold
  int thresh2=((struct thread_data*)arg)->thresh2;          //data threshold
  double fwhmy=((struct thread_data*)arg)->fwhmy;          //fwhm
  int boxCent=((struct thread_data*)arg)->boxCent;      //box size for mpfit
  int nmax=((struct thread_data*)arg)->nmax;          //data threshold
  int maxIt=((struct thread_data*)arg)->maxIt;          //data threshold

  struct cand_point *cand_head=((struct thread_data*)arg)->cand_list;  //output list of found points
  int verbose=((struct thread_data*)arg)->verbose;

  int boxsize;
  struct cand_point *cand_list;
  
  if(verbose == 1)
    {
      printf("In Thread\n"); 
      printf("Starting Processing\n"); 
    }
  /*---------------------------------------------------------------*/

  if(verbose == 1)
  
    {
    printf("Finding Points\n"); 
    }

  
  int np;
  //printf("%d %d %d %d %d %d %d\n",thresh1,thresh2,n_x,n_y,boxFind,nmin,nmax);

  cand_list=getRegions(image,thresh1,thresh2,boxFind,boxCent,n_x,n_y,nmin,nmax,imagemask,&np,verbose);

  /* get the automatic paramters if desired */
  if((fwhmx==0) || (fwhmy==0))
    {
      //printf("G %d %d %d %d %d %d %d\n",thresh1,thresh2,n_x,n_y,boxFind,nmin,nmax);

      getParams(cand_list,&fwhmx,&fwhmy);
    }

  //printf("XX %lf %lf\n",fwhmx,fwhmy);
  
  /*---------------------------------------------------------------*/
  
  if(verbose == 1)
    {
      printf("Found %i points\n",np); 
      printf("Starting Centroiding \n"); 
    }

  /*---------------------------------------------------------------*/

  /*Now cycle through the guesses. first do round/sharp checks a la daofind,
    if the point passes, then do the centroiding and add it to the list*/

  

  /*cycle through the possible detections. idl/fortran codes in
  original daofind uses lots of goto replaced with if statements*/
  
  cand_curr=cand_list;
  double *centroidVal;

  int iii=0;
  while(cand_curr!=NULL)
  {
    //get the center points of the possible detection, and its value
  
    centroidVal=windowedPos(image,cand_curr->x,cand_curr->y, boxCent,fwhmx,fwhmy,maxIt,n_x,n_y,verbose);
    cand_curr->x=centroidVal[0];
    cand_curr->y=centroidVal[1];
    cand_curr=cand_curr->next;
    free(centroidVal);
    iii=iii+1;
   	  //-----------------------------------------------------------//
	  
  } 

  if(verbose == 1)
    {
      printf("Finished Filtering/Centroiding\n");
      printf("Tidying Results\n");
    }

  //now assign the list back to the thread variable to pass it back
  ((struct thread_data*)arg)->cand_list=cand_list;

  //free memory
  //free(image);
  free(imagemask);

  //exit the thread properly
  pthread_exit(0);
}



  /*-------------------------------------------------------------------------------------*/

struct centroids *centroid(int *image, int n_x, int n_y, int thresh1, int thresh2, double fwhmx, double fwhmy,int boxFind, int boxCent,int *np, int nmin, int nmax,int maxIt, int verbose)
{

  /*main routine. parses input values, calls setup routine, divides up the threads, 
    runs the threads, combines and filters the outputs */

  struct centroids *output; 

  int ii,jj,i,j;          //counter variabls
  //int fname_check=0,hmin_check=0,fwhm_check=0;  //flags to check input

  if(verbose == 1)
    {
      printf("Starting Program\n");
      printf("  n_x = %d\n",n_x);
      printf("  n_y = %d\n",n_y);
      printf("  boxFind = %d\n",boxFind);
      printf("  boxCent = %d\n",boxCent);
      printf("  fwhmx = %lf\n",fwhmx);
      printf("  fwhmy = %lf\n",fwhmy);
      printf("  thresh2 = %d\n",thresh1);
      printf("  thresh2 = %d\n",thresh2);
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
  struct cand_point *real_top;  //List f points

  double x1,y1;  //positions of first value
  double x2,y2;  //positions of first value
  double p1,p2;   //fluxes of first and second avalue
  double rs;     //square o fdistance between x1,y1 and x2,y2
  
  double rmin=boxFind+1;  //radius to look for duplicate/false points
  int deadfirst=0;    //flag that first node has been deleted, to keep track fo pointers
  int firstval=1;     //flag that we're looking at the first node, to keep track of pointers
  char filename[sizeof "file100.fits"];

  /*------------------------------------------------------------*/
  
  int boxsize=maxValI(boxCent,boxFind);

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

	  thread_data_array[ind].thresh1=thresh1;        
	  thread_data_array[ind].thresh2=thresh2;        
	  thread_data_array[ind].fwhmx=fwhmx;        
	  thread_data_array[ind].fwhmy=fwhmy;        
	  thread_data_array[ind].boxFind=boxFind;
	  thread_data_array[ind].boxCent=boxCent;
	  thread_data_array[ind].nmin=nmin;
	  thread_data_array[ind].nmax=nmax;
	  thread_data_array[ind].maxIt=maxIt;
	  thread_data_array[ind].cand_list=cand_list[ind];
	  thread_data_array[ind].verbose=verbose;

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
    int iii,jjj;
    int inloop=0;
    cand_val=NULL;
    iii=0;

    for(ii=0;ii<NTHREAD-1;ii++)
      {

	iii=0;
	jjj=0;
	int skipit=0;
	if(ii > 0 & inloop ==1)
	  {
	    if(thread_data_array[ii].cand_list == NULL)
	      {
		skipit=1;
	      }
	    else
	      {
		//and get the next list
		cand_val->next=thread_data_array[ii].cand_list;
		cand_val=cand_val->next;
	      }
	  }


	if(inloop==0)
	  {
	    cand_val=thread_data_array[ii].cand_list;
	  }

	if(skipit==0)
	  {
	if(cand_val != NULL)  //check for empty list
	  {

	    //mark the start of the list (incase first part is null)
	    if (inloop==0)
	      {
		top_val=cand_val;
		inloop=1;
	      }
	     
	      
	    while(cand_val->next != NULL)
	      {

		if(verbose==1)
		  {
		    //printf("BB %d %lf %lf\n",ii,cand_val->x,cand_val->y);
		  }

		//Add the offset for that subimage


		cand_val->x=cand_val->x+thread_data_array[ii].fpix0-1;
		cand_val->y=cand_val->y+thread_data_array[ii].fpix1-1;


		cand_val=cand_val->next;
		iii=iii+1;

	      }
	    //And the last point in the list
	
	    printf("DD %d %d %lf %lf %d\n",ii,iii,cand_val->x,cand_val->y,thread_data_array[ii].fpix1-1);

	    cand_val->x=cand_val->x+thread_data_array[ii].fpix0-1;
	    cand_val->y=cand_val->y+thread_data_array[ii].fpix1-1;
		
	    jjj=1;
	    if(fabs((cand_val->x-(int)cand_val->x) - 0.642373) < 0.00001)
	      {
		printf("LL %d %d %lf %lf %d %d\n",iii,jjj,cand_val->x,cand_val->y,thread_data_array[ii].fpix0-1,thread_data_array[ii].fpix1-1);
		
	      }
	    if(fabs((cand_val->x-(int)cand_val->x) - 0.848779) < 0.00001)
	      {
		printf("MM %d %d %lf %lf %d %d\n",iii,jjj,cand_val->x,cand_val->y,thread_data_array[ii].fpix0-1,thread_data_array[ii].fpix1-1);
		
	      }

	    
	      if(iii==1)
	      {
	      	cand_val->x=cand_val->x-thread_data_array[ii].fpix0-1;
	      	cand_val->y=cand_val->y-thread_data_array[ii].fpix1-1;
	      }
	      
	if(iii==1)
	  {
	    printf("here\n");
	  }
	printf("JJ %d\n",iii);
	  }
	  }
      }
    /*and the same for the last segment (or a single segment in the unthreaded case, when
      the above loop is not executed)*/

    
    if(NTHREAD==1)
      {
	cand_val=thread_data_array[0].cand_list;
	top_val=cand_val;
      }


    if(cand_val != NULL)  //Are there items in the list?
      {

	//and get the next list

	if(NTHREAD > 1)
	  {
	    cand_val->next=thread_data_array[ii].cand_list;
	    cand_val=cand_val->next;
	    
	    
	    
	    while(cand_val->next != NULL)
	      {
		if(verbose==1)
		  {
		    //printf("ZZ %lf %lf\n",cand_val->x,cand_val->y);
		  }
		cand_val->x=cand_val->x+thread_data_array[NTHREAD-1].fpix0-1;
		cand_val->y=cand_val->y+thread_data_array[NTHREAD-1].fpix1-1;
		iii=iii+1;
		cand_val=cand_val->next;
	      }
	  }

      }
    cand_val->x=cand_val->x+thread_data_array[ii].fpix0-1;
    cand_val->y=cand_val->y+thread_data_array[ii].fpix1-1;

    //Now filter out points with badly failed fits (fqual > 5)

    curr_val=top_val;
    curr_pre=top_val;

    while(curr_val != NULL)
    {

      if((curr_val->y < 1000) | (curr_val->x < 1000))
      {
	printf("!!! %lf %lf \n",curr_val->x,curr_val->y);
      }

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
      as there is no previous nosde*/
    
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
      output[i].x2=curr_val->x2;
      output[i].y2=curr_val->y2;
      output[i].peak=curr_val->peak;
      output[i].xy=curr_val->xy;
      output[i].qual=curr_val->qual;
      curr_val=curr_val->next;
    }

  freeAll(&top_val);
		     
		       
  
  np[0]=ii;

  for (ii=0;ii<XSPLIT;ii++)
    {
      for (jj=0;jj<YSPLIT;jj++)
	{

	ind=ii+XSPLIT*jj;
	free(thread_data_array[ind].image);
	}
    }
      
  return output;

}

/*--------------------------------------------------------*/
