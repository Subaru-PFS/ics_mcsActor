
/* routines called by centroid, put into a separate file to make
   version control easier for different versions.*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h> 
#include "mpfit.h"
#include "centroid_types.h"
#include "fitsio.h"
#include "centroid.h"
#include "fitsio.h"


void printresult(double *x, mp_result *result,int num) 
{
/* Quick routine to print the results for a particular mpfit - label,
   number of iterations and the final value of the parameters. for
   debugging only */

  int i;

  if ((x == 0) || (result == 0)) return;  //NO RESULTS

  printf("%d ",result->niter);
  for (i=0; i<result->npar; i++) 
    {
      printf(" %lf ",x[i]);
    }
  printf("\n");


}
//--------------------------------------------------//
int gaussfunc1d(int m, int n,double *p, double *dy, double **dvec, void *vars)
{
  /* 
   * 1D Gaussian function for mpfit
   *
   * m - total number of data points
   * n - total number of parameters
   * p - array of length n with the parameters
   * dy - array of residuals to be returned
   * dvec - optional direct calculation of partial derivatives (not used here)
   * vars - private data to be passed (coordinates, data, etc) (struct vars_struct *)
   *
   * RETURNS: error code (0 = success)
   */

  int i;    //counter
  struct vars_struct *v = (struct vars_struct *) vars;  //structure to pass private data
  int *flux;  //data values
  int *ferr;  //errors
  double xc;        //the value of the centre points of the image

  //set the local variables to the structure.

  flux = v->flux;
  ferr = v->ferr;

  //cycle through the values.

  for (i=0;i<m;i++)
    {
      
      xc=i-p[0];

      //equation 

      dy[i] = ((double)flux[i] - p[2]*exp(-0.5*(xc*xc/(p[1]*p[1]*0.180337)))-p[3])/(double)ferr[i];

    }
  return 0;
}
//--------------------------------------------------//
int gaussfunc2d(int m, int n, double *p, double *dy, double **dvec, void *vars)
{
  /* 
   * 2D Gaussian function for mpfit
   *
   * m - total number of data points
   * n - total number of parameters
   * p - array of length n with the parameters
   * dy - array of residuals to be returned
   * dvec - optional direct calculation of partial derivatives (not used here)
   * vars - private data to be passed (coordinates, data, etc) (struct vars_struct *)
   *
   * RETURNS: error code (0 = success)
   */

  int i,j,n1;  //counters, size of image
  struct vars_struct *v = (struct vars_struct *) vars;  //data type for private variables
  int *flux, *ferr;  //data and errors
  double xc, yc;        //value of the centroids

  //Set the local variables to the structure

  flux=v->flux;
  ferr = v->ferr;

  //assuming a square image, get the dimensions of each side

  n1=(int)sqrt(m);

  //Cycle through the values. The data/residuals are 1D,
  //Map the coordinates to that assumign a square

  for (i=0;i<n1;i++)
    {
      for (j=0;j<n1;j++)
	{
	  //centre values
	  xc=i-p[0];
	  yc=j-p[1];

	  //Equations assuming independent FWHM in X and Y directions

	  dy[i*n1+j] = (((double)flux[i*n1+j] - p[4]*exp(-0.5*(xc*xc/(p[2]*p[2]*0.180337)+yc*yc/(p[3]*p[3]*0.180337)))-p[5]))/ferr[i*n1+j];
	  
	}
    }

  return 0;
}
//--------------------------------------------------//

int *make_mask(int* image,  int nbox, long n_x, long n_y, double mthresh,int VERBOSE)
{

  /*routine to make a mask for the image to pass to the convolution
    routine. Marks points that are within nbox pixels of a point above
    the threshold */

  int ii,jj,i,j;  //counters
  int *imagemask = calloc(n_x*n_y,sizeof(int));;

  //cycle through the pixels

  for (i=nbox;i<n_y-nbox;i++)
    {
      for (j=nbox;j<n_x-nbox;j++)
	{

	  //is the point currently unset, and above the threshold
	  if(imagemask[i*n_x+j]==0 && image[i*n_x+j] > mthresh)

	  {
	  
          //expand the region slightly to make sure that any pixel that will
          //be used in the convolution or peak finding is flagged.
	    for(ii=-nbox;ii<nbox;ii++)
	      {
	  	for(jj=-nbox;jj<nbox;jj++)
	  	  {
	  	    imagemask[(i+ii)*n_x+(j+jj)]=1;
	  	  }
	      }
	  }
      }
    }
  
  return imagemask;

} 
/*------------XXX--------------------------------------------*/

int *convol_sep(int *image, int *imagemask, double *kernel, int n_x, int n_y, int nbox,int VERBOSE)
{

  /*take an image and do a simple convolution with a kernel


    only pixels that are flagged in imagemask are convolved, as the 
    image consists mostly of empty space

    assumes the kernel is square

    as the gaussian kernel is separable, this is split into two 1d
    convolutions.

    note that for a small kernel compared to the image, direct
    convolution is faster than FFT.

    image - input image
    kernel - input kernel
    outimate - output image
    n_x, n_y - size of image
    nbox - size of kernel (assuming square)
 */

  int i,j,k,kk; // counter variables
  int *outimage;  //output image
  int *outimage1;  //temporary image

  outimage = calloc(n_x*n_y,sizeof(int));
  outimage1 = calloc(n_x*n_y,sizeof(int));
  
  int kcent=nbox/2;  //half size of box

  //cycle through the first time
  for (j=nbox;j<n_y-nbox;j++)
    {
      for (i=nbox;i<n_x-nbox;i++)
      {  	
  	
	//if a valid point, do the convolution
  	if(imagemask[j*n_x+i]==1)
  	  {
  	    for (k=0;k<nbox;k++)
  	      {
  		kk = k - kcent;
  		outimage[j*n_x+i] += (int)(image[j*n_x+(i+kk)]*kernel[k]);
  	      }
  	  }
      }
    }
  
  //now the same in the other direction 
  for (i=nbox;i<n_y-nbox;i++)
    {
      for (j=nbox;j<n_x-nbox;j++)
      {
  
  	if(imagemask[i*n_x+j]==1)
  	  {
  	    for (k=0;k<nbox;k++)
  	      {
  		kk = k - kcent;
  		outimage[i*n_x+j] += (int)(image[(i+kk)*n_x+(j)]*kernel[k]);
  	      }
  	    }
      }
    }
  


  //free memory and return
  free(outimage1);
  return outimage;

  }
/*-----------------------------------------------------------*/
int input_values(int *nbox1,double **gx,double **c1,int **mask, int *pixels,int *nhalf,double fwhm,int VERBOSE)
{

  /*calculates the various parameters used in the finding poriton - the kernel,
    the scaled kernel, the pixel mask, the centre point and the number of pixels*/

  //get nbox
  int nbox=2*floor(0.637*fwhm) + 1;   //boxsize for smoothing

  *nbox1=nbox;

  //Allocate memory for kernel information
  //*gx = malloc(nbox*sizeof(double));  
  *gx = calloc(sizeof(double),nbox);
  *c1 = malloc(nbox*sizeof(double));  
  *mask = malloc(nbox*nbox*sizeof(int));

  //values for gaussian
  double sigsq = (fwhm/2.35482 )*(fwhm/2.35482 );  // square of FWHM in sigma terms
  double radsq=(0.637*fwhm)*(0.637*fwhm);  // square of the radius of the fitting area

  int i,j;   //counters

  double rval;    //radius from centre of the fitting region
  double x1,y1;  // variable for x/y pixel position
  double gxtot=0; //for normalizing the kernel (1d)

  double row2[nbox],sumc1=0,sumc1sq=0;
  (*nhalf)=floor(0.637*fwhm); // centre point of box region for find

  if(VERBOSE==1)
    {
      printf("Creating Convolution Kernel\n");
      printf("nbox=%i\n",nbox);
      printf("nhalf=%i\n",(*nhalf));
    }

  //calculate 1d, 2d kernels, and masks

  for (i = 0;i< nbox;i++)
    {
      x1=i-(*nhalf);                 
      (*gx)[i]=exp(-0.5*x1*x1/sigsq);   //the separated kernel
	  
      gxtot=gxtot+(*gx)[i];

      //mask out pixels at r > 1.5 sigma
     for (j=0;j<nbox;j++)
     	{
     
     	  y1=j-(*nhalf);
     
     	  rval = x1*x1+y1*y1;
     
     	  if(rval > radsq)
     	    {
     	      (*mask)[i*nbox+j] = 0;
     
     	    }
     	  else
     	    {
     	      (*mask)[i*nbox+j] = 1;
     	      (*pixels)=(*pixels)+1;
     
     	    }
     
     	}

    }

  //exclude the centre pixel and reduce the number of valid pixels by 1
  (*mask)[(*nhalf)*nbox+(*nhalf)]=0; 
  (*pixels) = (*pixels) -1;       

  //scaling the convolution kernel
  for (i = 0;i< nbox;i++)
    {
      (*gx)[i]=(*gx)[i]/gxtot;

    }
  //c1 is a scaled version of the convolution kernel, summed along one
  //axis. used for round
  
  for (i=0;i<nbox;i++)
    {
      row2[i]=(i-(*nhalf))*(i-(*nhalf));  //centre position squared
    }
  for (i = 0;i< nbox;i++)
    {
      (*c1)[i]=exp(-0.5*row2[i]/sigsq);     //gaussian value
      sumc1=sumc1+(*c1)[i];                 //sum of c1
      sumc1sq=sumc1sq+(*c1)[i]*(*c1)[i];       //sum of c1^2
    }

  //scale value
  sumc1=sumc1/nbox;                     
  sumc1sq=sumc1sq-sumc1;
      
  for (i = 0;i< nbox;i++)
    {
      (*c1)[i]=((*c1)[i]-sumc1)/sumc1sq;
    }

  return 0;
}
/*--------------------------------------------------------*/

double calculate_sharp(int *image,int *mask, int nhalf, int ix, int iy, int n_x, int n_y,int pixels,int nbox)
{

  /*Calculate the sharpness parameter*/

  int sharptot=0;  //temporary variable for calculating sharp
  int ii,jj;       //counter variables
  double sharp;    //sharp parameter

  //cycle through the box around the centre point, and calculate the
  //sharpness statistic.

  //ii and jj cycle through -nhalf to nhalf
  //ix +ii and iy +jj are the coordinates on the image
  //ii + nhalf and jj+nhalf are the coodinates within nbox
  
  for (ii=-nhalf;ii<=nhalf;ii++)
    {
      for (jj=-nhalf;jj<=nhalf;jj++)
	{
	  //multiplying by the mask
	  sharptot += (double)image[(iy+jj)*n_x+(ix+ii)]*mask[(jj+nhalf)*nbox+(ii+nhalf)]; 
	}
      
    }

  //now scale and calculate the statistic
  sharp=((double)image[iy*n_x+ix]-sharptot/pixels)/(double)image[iy*n_x+ix];

  return sharp;
}
/*--------------------------------------------------------*/

double calculate_round(int *image,int *mask, int nhalf, int ix, int iy, int n_x, int n_y,int pixels,int nbox,double *c1)
{

  /*Calculate the roundness parameter*/

  int ii,jj,i;     //counters
  double around;   //round paramter
  double dxt[nbox],dyt[nbox];  //variable for calculating round
  double dx=0,dy=0;            //variable for calculating round

  //initialize variables
  for (i=0;i<nbox;i++)
    {
      dxt[i]=0;
      dyt[i]=0;
    }

  //dxt is the image summed along an axis
  for (ii=-nhalf;ii<=nhalf;ii++)
    {
      for (jj=-nhalf;jj<=nhalf;jj++)
	{
	  //this is doing a sum in the columns and rows, respectively
	  //by swapping ii and jj in the image. 
	  dxt[ii+nhalf]=dxt[ii+nhalf]+(double)image[(iy+jj)*n_x+(ix+ii)];
	  dyt[ii+nhalf]=dyt[ii+nhalf]+(double)image[(iy+ii)*n_x+(ix+jj)];
	}
    }
  
	  
  //and now for a total - estimates of roundness in each direction
  
  //dx ix the sum of dxt times a scaled version of the kernle summed
  //along the same axis
  
  for (ii=0;ii<nbox;ii++)
    {
      dx=dx+dxt[ii]*c1[ii];
      dy=dy+dyt[ii]*c1[ii];
    }
  
  around = 2.*(dx-dy) / ( dx + dy );    //Roundness statistic
 
  return around;

}
/*--------------------------------------------------------*/

struct cand_point *getfit_simple(int *image,int hpoint, int n_x,int n_y, int ix, int iy,double fwhm,int boxsize)
{

  /*Simple weighted centroid of values above a threshold. Fastest, but
    least accurate*/

  //counter variables

  int ii,jj; //counter variables
  struct cand_point *cand_curr;  //holds the centroiding reusults to pass back

  double xc,yc; //centroid
  double ssum;  //sum of fluxes for weighting

  //threshold for the weighted centroid, to avoid 
  //including background noise in teh calculation.

  int vthresh=image[(iy)*n_x+(ix)]/3;

  //initialize variables
  xc=0;   
  yc=0;   
  ssum=0;

  cand_curr=(struct cand_point*)malloc(sizeof(struct cand_point));


  //cycle through a box around the centre point

  for (ii=-hpoint;ii<=hpoint;ii++)
    {
      for (jj=-hpoint;jj<=hpoint;jj++)
	{
	  if(image[(iy+jj)*n_x+(ix+ii)] > vthresh)
	    {
	      
	      //x centroid, y centroid, sum for weighting
	      xc=xc+ii*(double)image[(iy+jj)*n_x+(ix+ii)];
	      yc=yc+jj*(double)image[(iy+jj)*n_x+(ix+ii)];
	      ssum=ssum+(double)image[(iy+jj)*n_x+(ix+ii)];
	    }
	}
    }
  
  //normalize
  xc=xc/ssum;
  yc=yc/ssum;

  
  //assign the results. We have no information about FWHM/Peak/background
  //for this method. 

  cand_curr=(struct cand_point*)malloc(sizeof(struct cand_point));
  cand_curr->x=ix+yc;
  cand_curr->y=iy+xc;
  cand_curr->fx=0;
  cand_curr->fy=0;
  cand_curr->peak=0;
  cand_curr->bg=0;
  cand_curr->qual=0;

  return cand_curr;

}
/*--------------------------------------------------------*/

struct cand_point *getfit_2d(int *image,int hpoint, int n_x,int n_y, int ix, int iy,double fwhm,int boxsize)
{
  int ii,jj;  // couners

  struct cand_point *cand_curr;

  //mpfit variables

  struct vars_struct v; //private structure with data/function information
  mp_result result;     //structure with results
  double xc,yc; //centroid

  double perror[6];	//errors in returned parameters
  mp_par pars[6];	//var for mpfit, holds information about limiting/fixing parameters
  int *subreg;          //variables for passing to MPFIT
  int *ferr;
  double pp[6];      //array holding fitted values in MPFIT
  int status;
  int npoints=boxsize*boxsize;
  subreg = malloc(boxsize*boxsize*sizeof(int));  
  ferr = malloc(boxsize*boxsize*sizeof(int));    

  memset(&result,0,sizeof(result));  
  result.xerror = perror;
  memset(pars,0,sizeof(pars));    
  
  for (ii=-hpoint;ii<=hpoint;ii++)
    {
      for (jj=-hpoint;jj<=hpoint;jj++)
	{
	  
	  subreg[(jj+hpoint)*boxsize+(ii+hpoint)]=image[(iy+jj)*n_x+(ix+ii)];
	  ferr[(jj+hpoint)*boxsize+(ii+hpoint)]=1;
	}
    }
  v.flux = subreg;
  v.ferr = ferr;
  pp[0] = (double)hpoint;
  pp[1] = (double)hpoint;
  pp[2] = 3.;
  pp[3] = 3.;
  pp[4] = subreg[hpoint*boxsize+hpoint];
  pp[5] = subreg[0];

  status=mpfit(gaussfunc2d, npoints, 6, pp, pars, 0, (void *) &v, &result);

  xc=pp[1];
  yc=pp[0];

  cand_curr=(struct cand_point*)malloc(sizeof(struct cand_point));
  cand_curr->x=ix+xc-boxsize/2;
  cand_curr->y=iy+yc-boxsize/2;

  cand_curr->fx=pp[3];
  cand_curr->fy=pp[2];
  cand_curr->peak=pp[4];
  cand_curr->bg=pp[5];

  if(status != 1)
    {
      cand_curr->qual=100;
    }
  else if((cand_curr->fx < 2) || (cand_curr->fy < 2))
    {
      cand_curr->qual=200;
    }
  else
    {
    cand_curr->qual=sqrt((cand_curr->x-ix)*(cand_curr->x-ix)+(cand_curr->y-iy)*(cand_curr->y-iy));
    }
  free(subreg);
  free(ferr);
  
  return cand_curr;
}
/*--------------------------------------------------------*/

struct cand_point *getfit_1d(int *image,int hpoint, int n_x,int n_y, int ix, int iy,double fwhm,int boxsize)
{

  /* Routine to calculate the centroid. This algorithm is based on
     fitting three 1-D fits in slices across the peak, and then interpolating
     the centre position (needed if the PSF is oriented at an angle) 

     This routine is faster than the 2D fit of getfit_2d, but is more 
     fragile to bad data or stray pixels, and doesn't explicitly fit the 
     peak of the PSF. 

     The FWHM in each direction is the mean of the three linear fits,
     the peak is the pixel value at the centre. The background is the 
     average value for all six fits.  
  */


  //counter variables

  int ii,jj,iii;

  struct cand_point *cand_curr;

  //mpfit variables

  struct vars_struct v; //private structure with data/function information
  int mpfit_status;          //exit status
  mp_result result;     //structure with results
  double xc,yc; //centroid

  double perror[4];	//errors in returned parameters
  mp_par pars[4];	//var for mpfit, holds information about limiting/fixing parameters
  int *subreg;          //variables for passing to MPFIT
  int *ferr;
  double pp[4];      //array holding fitted values in MPFIT

  int npoints=boxsize;

  subreg = malloc(boxsize*sizeof(int));  
  ferr = malloc(boxsize*sizeof(int));    

  //variables for the holding the centroids for each slice

  int ncut=3;       //number of slices to fit. Three has been optimal in testing.
  double ncutf=3.;
  double xfitbox[3];   
  double yfitbox[3];

  double fx=0;  //estimate of FWHM in the x direction
  double fy=0;  //estimate of FWHM in the y direction
  double bg=0;  //estimate of background

  //variables used for the interpolation of the fit

  double a1,b1,a2,b2;  //fits, as in y=ax+b
  double xsum,ysum,x2sum,xysum;  //sum of x, sum of y sum of x^2, sum of x*y
  int ival; //pixel value in fit 

  fx=0;
  fy=0;
  bg=0;
	      
  //extract the subregion
  for (iii=-1;iii<=1;iii++)
    {
      for (ii=-hpoint;ii<=hpoint;ii++)
	{
	  subreg[ii+hpoint]=image[(iy+iii)*n_x+(ix+ii)];
	  ferr[ii+hpoint]=1.;
	}

      //set the memory for mpfit
      memset(&result,0,sizeof(result));  
      result.xerror = perror;
      memset(pars,0,sizeof(pars));    
      
      //assign the data to the private array
      
      v.flux = subreg;
      v.ferr = ferr;

      //guess for the peak is the peak in x data
      pp[2]=subreg[hpoint];
      pp[0]=hpoint;
      pp[1]=fwhm;
      pp[3]=subreg[0];

      //call the fit

      mpfit_status= mpfit(gaussfunc1d, npoints, 4, pp, pars, 0, (void *) &v, &result);
      //keep track of the results
      xfitbox[iii+1]=pp[0];

      //estimate of the background and FWHM
      fx=fx+pp[1]/ncutf;
      bg=bg+pp[3]/(2*ncutf);  
      
    }

  //same in the other direction
  for (iii=-1;iii<=1;iii++)
    {
      for (jj=-hpoint;jj<=hpoint;jj++)
	{
	  subreg[jj+hpoint]=image[(iy+jj)*n_x+(ix+iii)];
	  ferr[jj+hpoint]=1;
	}
      
      memset(&result,0,sizeof(result));  
      result.xerror = perror;
      memset(pars,0,sizeof(pars));    
      
      //assign the data to the private array
      
      v.flux = subreg;
      v.ferr = ferr;
      pp[2]=subreg[hpoint];
      pp[0]=hpoint;
      pp[1]=fwhm;
      pp[3]=subreg[0];
      
      //call the fit

      mpfit_status= mpfit(gaussfunc1d, npoints, 4, pp, pars, 0, (void *) &v, &result);
      
      yfitbox[iii+1]=pp[0];
      fy=fy+pp[1]/ncutf;
      bg=bg+pp[3]/(ncutf*2);
    }

  //interpolate the fits for the true centre. This is necessary if the PSF is not
  //oriented in the XY direction

  //first, we fit a line using the three peak values in the x direction
  //this is a basic least squares fit

  xsum=0;
  ysum=0;
  x2sum=0;
  xysum=0;

  for (ii=0;ii<ncut;ii++)
    {
      //ival=ii+boxsize/2-1-ncut/2;
      ival=ii-ncut/2;

      ysum+=xfitbox[ii];
      xsum+=(double)ival;
      xysum+=xfitbox[ii]*(double)ival;
      x2sum+=xfitbox[ii]=xfitbox[ii];
    }
  
  a1=(ncut*xysum-xsum*ysum)/(ncut*x2sum-xsum*xsum);
  b1=(ysum-a1*xsum)/ncutf;

  //the same fit in the y direction

  xsum=0;
  ysum=0;
  x2sum=0;
  xysum=0;

  for (ii=0;ii<ncut;ii++)
    {
      //ival=ii+boxsize/2-1-ncut/2;
      ival=ii-ncut/2;

      ysum+=yfitbox[ii];
      xsum+=(double)ival;
      xysum+=yfitbox[ii]*(double)ival;
      x2sum+=yfitbox[ii]*yfitbox[ii];
    }
  
  a2=(ncut*xysum-xsum*ysum)/(ncut*x2sum-xsum*xsum);
  b2=(ysum-a2*xsum)/ncutf;

  //calculate the intercept of the two lines. 

  yc=(a2*b1+b2)/(1-a1*a2);
  xc=(a1*b2+b1)/(1-a1*a2);	  

  //assign the results to the structure and pass it back we need to
  //scale xc and yc back to pixel coordinates - add the centre point
  //of the fitting and the coordinates of the initial guess. 

  cand_curr=(struct cand_point*)malloc(sizeof(struct cand_point));
  cand_curr->x=xc-boxsize/2+ix;
  cand_curr->y=yc-boxsize/2+iy;
  cand_curr->peak=subreg[hpoint];
  cand_curr->fx=fy;
  cand_curr->fy=fx;
  cand_curr->bg=bg;
  cand_curr->qual=2;

  free(subreg);
  free(ferr);

  return cand_curr;

}

/*--------------------------------------------------------*/
