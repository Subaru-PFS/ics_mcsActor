
3 python level functions defined in centroid_wrapper.pyx

get_homes_call 

  takes an image and returns the home positions

centroid_coarse_call
   
   takes an image, an arc image, and the home positions from the above,
   and returns centroids for a large move and returns the centroids
   matched with the homes

centroid_fine_call

   takes an image, the home positions, and the centroids from
   centroid_coarse_call as teh previous positions, and returns
   the centroids matched with the homes

The file test.py runs the above from home positions to large move to
small move. 

to compile the code

> python setup.py build_ext --inplace


to perform unit tests on the centroid routines

make -f Makefile.centroidtest
./test_centroid

For the moment, the director TestData, used for running
test.py and the unit tests is located at

https://drive.google.com/drive/folders/1iJ2wqxSEAKhi91snNrgOO5uA9fxp6kEI?usp=sharing
