
/* Script to run centroid code on an input image for memory leak checking and profiling */

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


  //image variables
  int n_x,n_y;  //image dimensions
  int *image; //the image
  char *datafile="test.fits";
  centroids *outval; //array of structures holding the results
  int np[1];

  /*----------------------------------------------------*/

  //test image read

  //read an image from file
  image=read_fits_int(datafile,&n_x,&n_y);
  struct centroids *output;
  output = centroid(image, n_x, n_y, 555,555,4.,9.,14,14,np,10,20,1);
 
  free(output);
  free(image);
  fscanf(stdin, "c"); 
  
}

