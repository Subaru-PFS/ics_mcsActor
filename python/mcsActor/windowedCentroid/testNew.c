#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h> 
#include "fitsio.h"
#include "centroid.h"
#include "centroid_types.h"
#include "fitsio.h"
#include "fitsio2.h"

int main()

{

  int n_x,n_y,i;
  int *image;
  
  image=read_fits_int("testFile.fits",&n_x,&n_y);

  struct cand_point *result;
  struct cand_point *result_head=NULL;
  
  int thresh1=2200;
  int thresh2=1200;
  int boxsize=6;
  int xsize=n_x;
  int ysize=n_y;
  int nmin=10;
  int nmax=90;
  int verbose=0;
  int VERBOSE=0;
  int *mask=calloc(sizeof(int),n_x*n_y);
  int npoints;

  result=getRegions(image,thresh1,thresh2,boxsize,xsize,ysize,nmin,nmax,mask,verbose,&np,VERBOSE);

  result_head=result;

  double *centroid;
  
  for (i=0;i<npoints;i++)
    {
      //printf("circle point %lf %lf\n",result_head->x,result_head->y);
      centroid=windowedPos(image,result_head->x,result_head->y,6,4,4,20,xsize,ysize,VERBOSE);
      //printf("%lf %lf %lf\n",centroid[0],centroid[1],centroid[2]);
      result_head=result_head->next;
      
    }
  return 0;


  
}
