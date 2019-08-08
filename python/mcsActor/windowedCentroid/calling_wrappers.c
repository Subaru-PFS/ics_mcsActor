
/* some wrapper routines to call the centroiding and fibre identification for closed loop tests */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h> 
#include "mpfit.h"
#include "centroid_types.h"
#include "centroid.h"
#include "fitsio.h"
#include "fibreid.h"



fibreid *centroid_coarse(int *image, int *arc_image,double *homes,int n_x,int n_y,int thresh1, int thresh2, double fwhmx, double fwhmy, int boxFind, int boxCent, int nhome, int nmin, int nmax, int maxIt int verbose)
{ 

  /* call the centroiding and finding for a large move*/

  double maxarm=58.; //for creating patrol region

  int i; //counter

  //setup the homepos array
  fibreid *homepos=malloc(sizeof(fibreid)*nhome);
  for (i=0;i<nhome;i++)
    {
      homepos[i].x=homes[i];
      homepos[i].y=homes[nhome+i];
      homepos[i].idnum=i+1;
    }

  //now the centroiding of the current position

  int np1[1];    //number of points detected
  int npval;     //size of spiral

  //call the centroiding
  centroids *positions=centroid(image,n_x, n_y,thresh1,thresh2,fwhmx,fwhmy,boxFind,boxCent,np1,nmin,nmax,maxIt,verbose);

  //assign the values to the currentpos array

  fibreid *currentpos=malloc(sizeof(fibreid)*np1[0]);
  for (i=0;i<np1[0];i++)
    {
      currentpos[i].x=positions[i].x;
      currentpos[i].y=positions[i].y;
      currentpos[i].idnum=i+1;
    }

  //create the data needed for the identification (should be moved to a higher level later)

  long *patrolregion=fibremap(homepos,n_x,n_y,maxarm,nhome);
  int *spiral=create_spiral(5,10,&npval);    

  int *adjacent_fibre_list=create_adjacent_fibre_list(homepos,nhome,maxarm);

  //and the fibre identification

  fibre_identify_coarse(homepos,currentpos,patrolregion,arc_image,spiral,adjacent_fibre_list,npval,n_x,n_y,np1[0],nhome);

  return homepos;


}

fibreid *centroid_fine(int *image, double *homes,double *xp,double *yp,int n_x,int n_y,int thresh1, int thresh2, double fwhmx, double fwhmy, int boxFind, int boxCent, int nhome, int nmin, int nmax, int maxIt int verbose)
{ 

  /* call the centroiding and finding for a small move. Here we need the previous positions from the previous iteration.*/

  double maxarm=60.; //for creating patrol region

  int i; //counter

  //setup the homepos array. 
  fibreid *homepos=malloc(sizeof(fibreid)*nhome);
  for (i=0;i<nhome;i++)
    {
      homepos[i].x=homes[i];
      homepos[i].y=homes[nhome+i];
      homepos[i].idnum=i+1;
      homepos[i].xp=xp[i];
      homepos[i].yp=yp[i];
    }

  //now the centroiding of the current position
  int np1[1];

  centroids *positions=centroid(image,n_x, n_y,thresh1,thresh2,fwhmx,fwhmy,boxFind,boxCent,np1,nmin,nmax,maxIt,verbose);

  fibreid *currentpos=malloc(sizeof(fibreid)*np1[0]);
  for (i=0;i<np1[0];i++)
    {
      currentpos[i].x=positions[i].x;
      currentpos[i].y=positions[i].y;
      currentpos[i].idnum=i+1;
    }

  //and the fibre identification

  long *patrolregion=fibremap(homepos,n_x,n_y,maxarm,nhome);

  fibre_identify_fine(homepos,currentpos,patrolregion,n_x,n_y,np1[0],nhome);

  return homepos;

}


