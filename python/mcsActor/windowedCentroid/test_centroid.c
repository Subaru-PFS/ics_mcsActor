
/* Unit and Integration tests for centroiding code. Tests each of the subroutines
   separately, then tests of the function call for the three different fit types.
   The file "testdata.fits" is required, and will be given a simple checksum
   test for. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "centroid_types.h"
#include "centroid.h"
#include "fitsio.h"
#include "fitswrappers.h"

int main(int argc, char **argv)

{

  /*----------------------------------------------------*/

  printf("Starting Unit Testing\n");

  //parameters for fitting

  int threshold=500;
  int boxsize=9;
  double fwhm=3;
  int fittype;  //type of fit

  //image variables
  int n_x,n_y;  //image dimensions
  int *image; //the image
  char *datafile="TestData/testdata.fits"; //Input file for testing

  //Expected Image Values

  int real_totimage=97158;
  int real_n_x=19;
  int real_n_y=15;

  //variables to contain the results 

  centroids *outval; //array of structures holding the results
  int np; //number of sources detected

  //counters
  int i,j; //counters
  //long n_x,n_y; //size of image

  //set verbose=1 for extensive debugging info to stdout
  int verbose=0;

  /*----------------------------------------------------*/

  //test image read

  //read an image from file
  image=read_fits_int(datafile,&n_x,&n_y);

  //simple checksum
  int totimage=0;
  for (i=0;i<n_x*n_y;i++)
    {
      totimage+=image[i];
  
    }

  //check dimensions
  if((real_n_x != n_x) || (real_n_y != n_y))
    {
      printf("  Image Dimensions are incorrect\n");
      printf("  Should be (%d,%d) are (%d,%d)\n",real_n_x,real_n_y,n_x,n_y);
    }

  //check image checksum
  else if (totimage != real_totimage)
    {
      printf("  Image sum is incorrect.\n");
      printf("  Should be %d is %d\n",real_totimage,totimage);
    }
  else
    {
      printf("  Test image read correctly.\n");
    }

  /*----------------------------------------------------*/

  //test variable initialization

  //variables used by input values, and expected values.

  int nbox;
  double *gx;
  double *c1;
  int *mask;
  int pixels;
  int nhalf;
  fwhm=3.;

  int real_nbox=3;
  int real_nhalf=1;
  int real_pixels=8;
  double real_gx[3]={0.297549,0.404902,0.297549};
  double real_c1[3]={-0.070319,0.140637,-0.070319};
  int real_mask[9]={1,1,1,1,0,1,1,1,1};
  double gtol=0.0001;
  int gxcheck=1;
  int c1check=1;
  int maskcheck=1;

  //calculate
  input_values(&nbox,&gx,&c1,&mask,&pixels,&nhalf,fwhm,verbose);

  //check values
  if(real_nbox != nbox)
    {
      printf("  Value of nbox incorrect.\n");
      printf("  Should be %d is %d\n",real_nbox,nbox);
    }
  else if(real_nhalf != nhalf)
    {
      printf("  Value of nhalf incorrect.\n");
      printf("  Should be %d is %d\n",real_nhalf,nhalf);
    }
  else if(real_pixels != pixels)
    {
      printf("  Value of pixels incorrect.\n");
      printf("  Should be %d is %d\n",real_pixels,pixels);
    }
  else
    {
      for (i=0;i<nbox;i++)
	{
	  gxcheck=gxcheck*(fabs(gx[i]-real_gx[i]) < gtol);
	  c1check=c1check*(fabs(c1[i]-real_c1[i]) < gtol);
	}
      for (i=0;i<nbox*nbox;i++)
	{
	  maskcheck=maskcheck*(mask[i]==real_mask[i]);
	}
     if(gxcheck==0)
	{
	  printf("  Variable gx is incorrect\n");
	  printf("  Should be ");
	    for (i=0;i<nbox;i++)
	      {
		printf("  %lf ",real_gx[i]);
	      }
	    printf("  \nis ");
	    for (i=0;i<nbox;i++)
	      {
		printf("  %lf ",gx[i]);
	      }
	}
     else if(c1check==0)
	{
	  printf("  Variable c1 is incorrect\n");
	  printf("  Should be ");
	    for (i=0;i<nbox;i++)
	      {
		printf("  %lf ",real_c1[i]);
	      }
	    printf("  \nis ");
	    for (i=0;i<nbox;i++)
	      {
		printf("  %lf ",c1[i]);
	      }
	}
     else if(maskcheck==0)
	{
	  printf("  Variable mask is incorrect\n");
	  printf("  Should be ");
	    for (i=0;i<nbox*nbox;i++)
	      {
		printf("  %d ",real_mask[i]);
	      }
	    printf("  \nis ");
	    for (i=0;i<nbox*nbox;i++)
	      {
		printf("  %d ",mask[i]);
	      }
	}

     else
       {
	 printf("  Results of input_values are correct\n");
       }

    }

  
  /*----------------------------------------------------*/

  //test make_mask wiht a simple checksum. 

  nbox=3;
  int hmin=500;
  int real_checkval=8694;

  int *imagemask=make_mask(image,nbox,n_x,n_y,hmin,verbose);
  int checkval=0;
  for (i=0;i<n_x*n_y;i++)
    {
      checkval+=i*imagemask[i];
    }
  
  if(checkval != real_checkval)
    {
      printf("  Checksum for imagemask is incorrect.\n");
    }
  else
    {
      printf("  Checksum for imagemask is correct.\n");
    }

  
  /*----------------------------------------------------*/


  //test convol_sep with a simple checksum

  nbox=3;
  hmin=500;
  real_checkval=8116017;
  int *h=convol_sep(image,imagemask,gx,n_x,n_y,nbox,verbose);
  checkval=0;
  for (i=0;i<n_x*n_y;i++)
    {
      checkval+=i*h[i];
      //printf("  %d ",imagemask[i]);
    }
  if(checkval != real_checkval)
    {
      printf("  Checksum for convolved image is incorrect.\n");
    }
  else
    {
      printf("  Checksum for convolved image is correct.\n");
    }

  /*----------------------------------------------------*/
  //check calculate_sharp and calculate_round

  int ix=9;
  int iy=7;
  nhalf=1;
  pixels=8;
  nbox=3;

  double real_sharp=0.324513;
  double real_around=0.048995;

  double sharp=calculate_sharp(image,mask,nhalf,ix,iy,n_x,n_y,pixels,nbox);
  double around=calculate_round(image,mask,nhalf,ix,iy,n_x,n_y,pixels,nbox,c1);

  if(fabs(sharp-real_sharp) > gtol)
    {
      printf("  calculate_sharp returns wrong value.\n");
      printf("  is %lf should be %lf\n",sharp,real_sharp);
    }
  else if(fabs(around-real_around) > gtol)
    {
      printf("  calculate_round returns wrong value.\n");
      printf("  is %lf should be %lf\n",around,real_around);
    }
  else
    {
      printf("  calculate_sharp and calculate_round return correct values\n");
    }

  /*----------------------------------------------------*/

  //check the fits, each called directly

  struct cand_point *cand_curr;
  struct cand_point *cand_curr1;
  int hpoint=(boxsize-1)/2;

  //fittype=0

  cand_curr=getfit_simple(image,hpoint,n_x,n_y,ix,iy,fwhm,boxsize);

  double x0=8.212636;
  double y0=7.407807;

  if((fabs(cand_curr->x-x0) > gtol) || (fabs(cand_curr->y-y0) > gtol))
    {
      printf("  centroid for fittype 0 is incorrect.\n");
      printf("  is (%lf,%lf) should be (%lf,%lf) \n",cand_curr->x,cand_curr->y,x0,y0);
    }
  else
    {
      printf("  centroid for fittype 0 is correct\n");
    }


  //fittype=1

  cand_curr=getfit_1d(image,hpoint,n_x,n_y,ix,iy,fwhm,boxsize);

  double x1=9.321869;
  double y1=6.618853;
  double fx1=2.834477;
  double fy1=2.749043;
  double peak1=1290.000000;
  double back1=346.013957;

    if((fabs(cand_curr->x-x1) > gtol) || (fabs(cand_curr->y-y1) > gtol))
    {
      printf("  centroid for fittype 1 is incorrect.\n");
      printf("  is (%lf,%lf) should be (%lf,%lf) \n",cand_curr->x,cand_curr->y,x1,y1);
    }
  else
    {
      printf("  centroid for fittype 1 is correct\n");
    }

    if((fabs(cand_curr->fx-fx1) > gtol) || (fabs(cand_curr->fy-fy1) > gtol) || (fabs(cand_curr->peak-peak1) > gtol) || (fabs(cand_curr->bg-back1) > gtol))
      {
	printf("  Parameters for fittype 1 are incorrect. For FWHMX, FWHMY, BG, Peak\n");
	printf("  is (%lf,%lf,%lf,%lf) should be (%lf,%lf,%lf,%lf) \n",cand_curr->fx,cand_curr->peak,cand_curr->bg,cand_curr->fy,fx1,fy1,peak1,back1);
      }
    else
      {
	printf("  Parameters for fittype 1 are correct\n");
      }


  //fittype=2
  double x2=9.406028;
  double y2=6.564925;
  double fx2=2.885531;
  double fy2=2.972002;
  double peak2=1100.932782;
  double back2=334.232987;



  cand_curr=getfit_2d(image,hpoint,n_x,n_y,ix,iy,fwhm,boxsize);

    if((fabs(cand_curr->x-x2) > gtol) || (fabs(cand_curr->y-y2) > gtol))
    {
      printf("  centroid for fittype 2 is incorrect.\n");
      printf("  is (%lf,%lf) should be (%lf,%lf) \n",cand_curr->x,cand_curr->y,x2,y2);
    }
  else
    {
      printf("  centroid for fittype 2 is correct\n");
    }

    if((fabs(cand_curr->fx-fx2) > gtol) || (fabs(cand_curr->fy-fy2) > gtol) || (fabs(cand_curr->peak-peak2) > gtol) || (fabs(cand_curr->bg-back2) > gtol))
      {
	printf("  Parameters for fittype 2 are incorrect. For FWHMX, FWHMY, BG, Peak\n");
	printf("  is (%lf,%lf,%lf,%lf) should be (%lf,%lf,%lf,%lf) \n",cand_curr->fx,cand_curr->fy,cand_curr->peak,cand_curr->bg,fx2,fy2,peak2,back2);
      }
    else
      {
	printf("  Parameters for fittype 2 are correct\n");
      }

    printf("Finished Unit Testing\n");
  printf("Starting Integration Testing\n");

  //call the routine from the top level, for each of the three fit types

  fittype=0;
  threshold=500.;
  fwhm=3.;
  verbose=0;
  boxsize=9;
  outval=centroid(image,n_x,n_y,threshold,fwhm,boxsize,&np,verbose,fittype);


  printf("  Testing fittype=0\n");
  if(np != 1)
    {
      printf("  Found the wrong number of sources.\n");
      printf("  Should be 1, is %d\n",np);
    }
  else 
    {
      checkval=1;
      checkval*=fabs(outval[0].x-x0) < gtol;
      checkval*=fabs(outval[0].y-y0) < gtol;

      if(checkval==1)
	{
	  printf("  centroid executed successfully for fittype=0\n");
	}
      else
	{
	  printf("  centroid  failed for fittype=0\n");
	  printf("  Should be (%lf,%lf), is  (%lf,%lf)\n",x0,y0,outval[0].x,outval[0].y);
	}
    }

  fittype=1;
  threshold=500.;
  fwhm=3.;
  verbose=0;
  boxsize=9;
  image=read_fits_int(datafile,&n_x,&n_y);
  outval=centroid(image,n_x,n_y,threshold,fwhm,boxsize,&np,verbose,fittype);


  printf("  Testing fittype=1\n");
  if(np != 1)
    {
      printf("  Found the wrong number of sources.\n");
      printf("  Should be 1, is %d\n",np);
    }
  else 
    {
      checkval=1;
      checkval*=fabs(outval[0].x-x1) < gtol;
      checkval*=fabs(outval[0].y-y1) < gtol;

      if(checkval==1)
	{
	  printf("  centroid executed successfully for fittype=1\n");
	}
      else
	{
	  printf("  centroid  failed for fittype=1\n");
	  printf("  Should be (%lf,%lf), is  (%lf,%lf)\n",x1,y1,outval[0].x,outval[0].y);
	}
    }


  fittype=2;
  threshold=500.;
  fwhm=3.;
  verbose=0;
  boxsize=9;
  image=read_fits_int(datafile,&n_x,&n_y);
  outval=centroid(image,n_x,n_y,threshold,fwhm,boxsize,&np,verbose,fittype);


  printf("  Testing fittype=1\n");
  if(np != 1)
    {
      printf("  Found the wrong number of sources.\n");
      printf("  Should be 1, is %d\n",np);
    }
  else 
    {
      checkval=1;
      checkval*=fabs(outval[0].x-x2) < gtol;
      checkval*=fabs(outval[0].y-y2) < gtol;

      if(checkval==1)
	{
	  printf("  centroid executed successfully for fittype=2\n");
	}
      else
	{
	  printf("  centroid  failed for fittype=1\n");
	  printf("  Should be (%lf,%lf), is  (%lf,%lf)\n",x2,y2,outval[0].x,outval[0].y);
	}
    }


  printf("Finished Integration Testing\n");

  return 0;
}



