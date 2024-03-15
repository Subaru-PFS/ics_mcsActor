#ifndef centroid_types_H_   /* Include guard */
#define centroid_types_H_

typedef struct QNode
{
  long key;
    struct QNode* next;
}QNode;

typedef 
struct Queue {
    struct QNode *front, *rear;
}Queue;

typedef struct thread_data
{
  //STRUCTURE TO HOLD THE INFORMATION TO BE PASSED TO AND FROM THE THREAD


  int *image;
  int n_x;
  int n_y;
  int fpix0;
  int fpix1;
  double thresh1;     //IMAGE THRESHOLD 
  double thresh2;     //IMAGE THRESHOLD 
  double fwhmx;     //FWHM
  double fwhmy;     //FWHM
  int boxFind;     //FWHM
  int boxCent;     //FWHM
  int nmin;
  int maxIt;
  int verbose;
  int np;

  struct cand_point *cand_list;  //LIST OF CANDIDATE POINTS (PASSED BACK TO ROUTINE)

} thread_data;


//--------------------------------------------------//

/* STRUCTURE TO HOLD A LIST OF INITIAL GUESSES (LOCAL MAXIMA ABOVE THRESHOLD)*/
typedef struct guess_point
{
  int x;   //X POSITION
  int y;   //Y POSITION
  struct guess_point *next;  
} guess_point;

//--------------------------------------------------//
typedef struct cand_point
{
  //RESULTS FROM MPFIT

  double x;  //X POSITION
  double y;  //Y POSITION
  double x2;  //FWHM X
  double y2;  //FWHM Y
  double peak;  //PEAK VALUE
  double xy; //ANGLE OF FWHM
  double qual; //quality flag
  double back; //number of points in region
  double npt; //number of points in region
  struct cand_point *next;
} cand_point;
//--------------------------------------------------//

typedef struct centroids
{
  //int npoints;
  double x;
  double y;
  double x2;
  double y2;
  double peak;
  double xy;
  double back;

} centroids;

typedef struct ag_point
{
  //RESULTS FROM AG ROUTINE

  double x;  //X POSITION
  double y;  //Y POSITION
  double fx;  //X POSITION
  double fy;  //Y POSITION
  double peak;  //PEAK VALUE
  double qual;
  struct ag_point *next;
} ag_point;

typedef struct ag_stars
{
  int npoints;
  double *x;
  double *y;
  double *fx;
  double *fy;
  double *peak;
  int *qual;
  int np;
} ag_stars;


/*MPFIT STRUCTURE - DEFINE THE PRIVATE STRUCTURE FOR THE DATA, ERRORS,
  COORDINATES ETC. ANY 2D BEHAVIOUR IS HERE.  REMOVED X AND Y BECAUSE WE
  CAN GET THESE, ASSUMING IMAGE IS SQUARE */

typedef struct vars_struct {
  int *flux;     //DATA
  int *ferr;     //ESTIMATE OF ERROR
} vars_struct ;

//--------------------------------------------------//
typedef struct fibreid
{
  double xp;        //list of previous x positions
  double yp;        //list of previous y position
  double xt;        //target x list
  double yt;        //target y list
  double xc;        //current x position
  double yc;        //current y position
  double x;        //centre coordinate
  double y;        //
  double peak;
  double back;
  double fx;
  double fy;
  int qual;
  int idnum;        //ID number

} fibreid;

//--------------------------------------------------//

//STRUCTURE TO HOLD A LIST OF INITIAL GUESSES (LOCAL MAXIMA ABOVE THRESHOLD)
typedef struct id_list
{
  //RESULTS FROM MPFIT
  int idnum;      //ID NUMBER
  int index;      //INDEX IN ORIGINAL STRUCTURE
  struct id_list *next;
  struct match_list *match;  //LIST OF POTENTIAL FIBRE MATCHES
} id_list;

typedef struct match_list
{
  //RESULTS FROM MPFIT
  int idnum;  
  struct match_list *next;
} match_list;



#endif // centroid_types_H_
