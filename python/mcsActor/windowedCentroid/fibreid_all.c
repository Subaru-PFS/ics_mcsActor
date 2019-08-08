
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "centroid_types.h"
#include "fitsio.h"
#include "fibreid.h"
#include <stddef.h>
    
void add_node(struct id_list **head,int val,int ind)
{

  /*add a node to an id_list*/

  //create node
  struct id_list *node;
  node = (struct id_list *) malloc( sizeof(struct id_list) ); 

  //update the values
  node->idnum=val;
  node->index=ind;
  node->match=NULL;

  //empty list case
  if(*head==NULL)
    {
      *head=node;
      node->next=NULL;
    }
  //add the node at the beginning
  else
    {
      node->next=*head;
      *head=node;
    }
}

void add_node_cur(struct id_list **head,int ind)
{
  /*create node*/

  struct id_list *node;
  node = (struct id_list *) malloc( sizeof(struct id_list) ); 

  //update the values
  node->index=ind;
  node->match=NULL;

  //empty list case
  if(*head==NULL)
    {
      *head=node;
      node->next=NULL;
      node->match=NULL;

    }
  //add at the beginning of the list
  else
    {
      node->next=*head;
      *head=node;
      node->match=NULL;
    }
}

void add_node_list(struct id_list **head,int val)
{
  /*create node*/

  struct match_list *node;
  node = (struct match_list *) malloc( sizeof(struct match_list) ); 
  struct match_list *node1;

  int dup=0;

  //update values
  node->idnum=val;

  //empty list case
  if((*head)->match==NULL)
    {

      (*head)->match=node;
      node->next=NULL;
    }
  //add at end of list
  else
    {

      //make sure it's not a duplicate
      node1=(*head)->match;

      while(node1 != NULL)
	{
	  if(node1->idnum==val)
	    {
	      dup=1;
	    }
	  node1=node1->next;
	}
      if(dup==0)
	{
	  node->next=(*head)->match;
	  (*head)->match=node;
	}
    }

}

void print_list_match(struct match_list *head)
{

  /*routine to print a match list*/

  struct match_list *ptr;       //temporary pointer

  //work through the list, printing each value
  ptr=head;
  while (ptr != NULL)
   {
     printf(" %d ",ptr->idnum);
     ptr=ptr->next;
   }
  printf("\n");
}

int print_list(struct id_list *head)
{

  //routine to print the contents of an ID list

  struct id_list *ptr;  //temporary pointer
  struct match_list *ptr1;

  int i=0;
  //work through the list
  ptr=head;
  while (ptr != NULL)
   {
     printf("%d %d\n",ptr->idnum,ptr->index);
     i=i+1;

     ptr1=ptr->match;

     print_list_match(ptr->match);

     ptr=ptr->next;
   }

  return i;
}


int remove_node(struct id_list **head,int val)
{

  /*routine to remove */

  struct id_list *temp,*prev;  //temporary pointers
  int val1;
  temp=*head;

  while (temp != NULL)
    {
      if(temp->idnum==val)
	{
	  val1=temp->index;

	  if(temp==*head)
	    {
	      *head=temp->next;
	    }
	  else
	    {
	      prev->next=temp->next;
	    }
	  free(temp);
	  return(val);
	    
	}
      else
	{
	  prev=temp;
	  temp=temp->next;
	}

    }

  return -1;
}
    

int remove_node_list(struct match_list **head,int val)
{

  /*remove node from a list*/

  struct match_list *temp,*prev;
  
  temp=*head;

  while (temp != NULL)
    {
      if(temp->idnum==val)
	{
	  if(temp==*head)
	    {
	      *head=temp->next;
	      free(temp);
	      return 1;
	    }
	  else
	    {
	      prev->next=temp->next;
	      free(temp);
	      return 1;
	    }
	}
      else
	{
	  prev=temp;
	  temp=temp->next;
	}

    }
  return 0;
}
    

int remove_node_cur(struct id_list **head,int ind)
{

  /*remove node from list*/

  struct id_list *temp,*prev;
  
  temp=*head;

  while (temp != NULL)
    {
      if(temp->index==ind)
	{
	  if(temp==*head)
	    {
	      *head=temp->next;
	      free(temp);
	      return 1;
	    }
	  else
	    {
	      prev->next=temp->next;
	      free(temp);
	      return 1;
	    }
	}
      else
	{
	  prev=temp;
	  temp=temp->next;
	}

    }
  return 0;
}

void xy_from_thetaphi(double l1,double l2,double theta,double phi,double xc,double yc,double *xf,double *yf)
{

  *xf=l1*cos(2*M_PI-theta)+l2*cos((2*M_PI-theta)+(M_PI-phi))+xc;
  *yf=l1*sin(2*M_PI-theta)+l2*sin((2*M_PI-theta)+(M_PI-phi))+yc;
}

long *fibremap(struct fibreid *homepos,int nx,int ny,double maxarm,int nhome)

/* 

   program to create a map of pixel positions which can be reached by
   more than one fibre

   INPUT

   homepos - home positions of fibres, in pixel coordinates
   nx,ny   - size of image

   OUTPUT

   an integer map of pixel position flags


 */

{

  long i,j,k; //counters

  long *patrolregion;   //map of patrol regions
  patrolregion=calloc(nx*ny,sizeof(long));

  double d_from_home;  //distance from home
  double maxarm_adj=maxarm+1;  //adjustved value of maximum arm length

  //for each home position
  for (i=0;i<nhome;i++)
    {

      //for each pixel within the potential maxarm region
      for (j=-maxarm_adj;j<=maxarm_adj;j++)
	{
	  for (k=-maxarm_adj;k<maxarm_adj;k++)
	    {
	      //calculate the distance in pixels from home
	      d_from_home=sqrt(j*j+k*k);

	      //if this is within the movement region for this fibre, increment. 
	      if(d_from_home <= maxarm_adj)
		{
		  //if it has no region yet, add ID_num
		  if(patrolregion[(j+(int)homepos[i].x)*ny+(k+(int)homepos[i].y)]==0)
		    {    
		      patrolregion[(j+(int)homepos[i].x)*ny+(k+(int)homepos[i].y)]+=homepos[i].idnum;
		    }
		  //if it already has two, add 100000000*idnum, to keep them separate

		  else if(patrolregion[(j+(int)homepos[i].x)*ny+(k+(int)homepos[i].y)] > 10000)
		    {
		      patrolregion[(j+(int)homepos[i].x)*ny+(k+(int)homepos[i].y)]+=homepos[i].idnum*100000000L;

		    }
		  //if it already has one, add 10000*idnum, to keep them separate
		  else
		    {
		      patrolregion[(j+(int)homepos[i].x)*ny+(k+(int)homepos[i].y)]+=homepos[i].idnum*10000L;

		    }
		}
	    }

	}
    }
  //and return
  return patrolregion;
}

/*--------------------------------------------------------------------*/
int is_hidden(int *arc_image,long nx,long ny,int xc,int yc,int inner_r, int outer_r,int thresh,float maxd)
{

  /* routine to check for hidden fibres. Looks around the position of
     the spots for pixels above the threshold, checks the separation
     between those pixels.
   */


  int i,j; //counters

  int npix=(int)M_PI*(outer_r*outer_r); //size of region
  int xpos[npix];  //x and y positions in region
  int ypos[npix];
  float dc; //distance from centre
  float dd; //distances betwee points
  int nmax=0;  //number of points above threshold

  //float maxd=5; //separation between bo
  int npairs; //number of separated pairs
  int hidden; //hidden flag

  //go through the box and find any maxima over 
  //the threshold, inside the radius outer radius

 
  for (i=-outer_r;i<=outer_r;i++)
    {
      for (j=-outer_r;j<=outer_r;j++)
	{

	  //calculate distance from centre

	  dc=sqrt(i*i+j*j);

	  //is the value > the threshold and in a circular radius?
	  //then add to the list and inside radius

	  if ((arc_image[(i+xc)*nx+(j+yc)] > thresh) && (dc < outer_r))
	    {

	      xpos[nmax]=i;
	      ypos[nmax]=j;
	      nmax+=1;
	    }  	           
	}
    }


  npairs=0;

  //no points found - the fibre isn't hidden
  if(nmax <= 2)
    {
      hidden=0;
    }

  //points found
  else 
    {
      hidden=1;
      //cycle through unique pairs of points , 
      for (i=0;i<nmax;i++)
	{
	  for (j=i+1;j<nmax;j++)
	    {
	      
	      //calculate distances
	      dd=sqrt((xpos[i]-xpos[j])*(xpos[i]-xpos[j])+(ypos[i]-ypos[j])*(ypos[i]-ypos[j]));
		//if the distance between points is > xposmu, the
		//fibre did not go behind the mask, or went behind and
		//came out again, and is not hidden
		if(dd > maxd)
		  {
		    hidden=0;
		    npairs+=1;
		  }
	    }	 	    
	}
    }
     

  return hidden;


}

/*--------------------------------------------------------------------*/

int *wherefrom(int *arc_image,long n_x,long n_y,int xc,int yc,int *spiral,int np)
{

  /*Routine to check what direction a fibre came from, from the arc image. Searches
   around the position of a fibre to find the arc, returns the position.*/

  int thresh=120;
  int i=0; //counter

  int notfound=1;  //flag for finding

  int *position;  //array containing position
  position=malloc(2*sizeof(int));

  
  int xx,yy;  //current x and y position for seraching
  //cycle through the search pattern
  while(notfound && i<np)
    {
      xx=xc+spiral[i*2];
      yy=yc+spiral[i*2+1];
      //is the pixel above the threshold
      if(arc_image[(xx)*n_x+(yy)] > thresh)
	{
	  notfound=0;
	  position[0]=xx;
	  position[1]=yy;
	}
      i=i+1;
    }

  //case where it's not found, set position to 0,0
  if(i==np)
    {
      position[0]=0;
      position[1]=0;

      return position;
    }
  else
    {
      return position;
    }
}


/*--------------------------------------------------------------------*/
int *create_spiral(int minr,int maxr,int *nval)
{
  /*generate a spiral pattern with a given inner and outer radius, for
    searching around a position*/

  int i; //counter
  int temp; //temporary variables

  int x=0,y=0,dx=0,dy=-1,np=0; //variables for generating the spiral

  int *spiral; //contain the spiral

  //the below calculates a spiral pattern around the point. First, get the size of the spiral.

  for (i=0;i<(2*maxr+1)*(2*maxr+1);i++)
    {
      if (((x*x+y*y)< maxr*maxr) && ((x*x+y*y)> minr*minr))
	{
	  np+=1;
	}
      if ((x == y) || ((x < 0) && (x == -y)) || ((x > 0) && (x == 1-y)))
	{
	  temp=dx;
	  dx=-dy;
	  dy=temp;

	}
      x=x+dx;
      y=y+dy;

    }
  
  //allocate the memory
  spiral=malloc(np*2*sizeof(int));

  x=0;y=0;dx=0;dy=-1;np=0;
  for (i=0;i<(2*maxr+1)*(2*maxr+1);i++)
    {
      if (((x*x+y*y)< maxr*maxr) && ((x*x+y*y)> minr*minr))
	{
	  
	  spiral[np*2]=x;
	  spiral[np*2+1]=y;


	  np+=1;
	}
      if ((x == y) || (x < 0 && x == -y) || (x > 0 && x == 1-y))
	{
	  temp=dx;
	  dx=-dy;
	  dy=temp;
	}
      x=x+dx;
      y=y+dy;

    }

  //set the variable with the size
  *(nval)=np;

  return spiral;

}

/*--------------------------------------------------------*/

int adjacent_assigned(struct id_list *un_fibre,long index,int *adjacent_fibre_list)
{

  /*routine to check if the fibres adjacent to the central fibre are assigned or not*/

  int i;  //counter
  long idnum; //id number fo the fibre (place in array)
  int idval; //id value for the fibre

  struct id_list *ptr_f; //pointer to ID list

  //check the six possible adjacent fibres
  for (i=0;i<6;i++)
    {
      //get the ID val
      idval=adjacent_fibre_list[index*6+i];

      //cycle through the unassigned fibres
      ptr_f=un_fibre;
      while(ptr_f != NULL)
	{

	  //IS the ID in the unassigned fibre list?
	  idnum=ptr_f->idnum;

	  //if so, return
	  if(idnum==idval)
	    {
	      return 0;
	     
	    }
	  else
	    {
	      ptr_f=ptr_f->next;
	    }

	}

    }

  //free the value
  free(ptr_f);

  //return if haven't found the fibre
  return 1;


}

/*--------------------------------------------------------*/

int *create_adjacent_fibre_list(struct fibreid *homepos,int nhome,double maxarm)
{

  /*
    create an array with list of adjacent COBRAs
 */
  float separation=maxarm*2;
  int *adjacent_fibres;
  
  int i,j,ii;
  float dd,xd,yd;

  //allocate memory. First element is the fibre id, 
  //next six are the adjacent ones

  adjacent_fibres=calloc(7*nhome,sizeof(int));

  //for each fibre
  for (i=0;i<nhome;i++)
    {
      ii=0;
      adjacent_fibres[i*7]=homepos[i].idnum;
	
      //look for centre positions close to it
      for (j=0;j<nhome;j++)
	{
	  xd=homepos[i].x-homepos[j].x;
	  yd=homepos[i].y-homepos[j].y;
	  
	  dd=sqrt((xd*xd+yd*yd));

	  //and assign
	  if((dd < separation) && (i!=j))
	    {
	      adjacent_fibres[i*7+ii+1]=homepos[j].idnum;
	      ii+=1;
	    }
	}
    }

  return adjacent_fibres;
}

/*--------------------------------------------------------------------*/

struct fibreid *read_homes(char *filename,int npos)
{

  /*routine to read COBRA centre positions from a file, for testing purposes.
   file is a text file with npos line, two coordinates per line.*/

  int i; //counter

  struct fibreid *homepos;  //structure for values
  
  //initialze memory

  homepos=malloc(sizeof(fibreid)*npos);

  double xx,yy; //positions

  FILE* file=fopen(filename,"r");

  for (i=0;i<npos;i++)
    {
      fscanf(file,"%lf %lf",&xx,&yy);
      homepos[i].y=xx;
      homepos[i].x=yy;
      homepos[i].idnum=i;
    }

  return homepos;

}
/*--------------------------------------------------------*/

void update_list(struct id_list **as_fibre,struct id_list **un_fibre,struct id_list **un_current,int idnum,int index,int ind,int indm,fibreid *homepos,fibreid *currentpos)
{

  //assign to the homepos structure. 

  homepos[index].xc=currentpos[indm].x;
  homepos[index].yc=currentpos[indm].y;

  //add the value to the list of assigned fibres
  add_node(as_fibre,idnum,index);
  //remove it from the list of unassigned fibres
  remove_node(un_fibre,idnum);
  //remove ind from list of unassigned positions
  remove_node_cur(un_current,ind);

  return; 
}
/*--------------------------------------------------------*/

void update_list_hidden(struct id_list **as_fibre,struct id_list **un_fibre,struct id_list **un_current,int idnum,int index,double xc,double yc,fibreid *homepos,fibreid *currentpos)
{

  //assign to the homepos structure. 

  homepos[index].xc=xc;
  homepos[index].yc=yc;

  //add the value to the list of assigned fibres
  add_node(as_fibre,idnum,index);
  //remove it from the list of unassigned fibres
  remove_node(un_fibre,idnum);
  //remove ind from list of unassigned positions

  return; 
}
/*--------------------------------------------------------*/

void  initial_position(struct id_list **as_fibre,struct id_list **un_fibre,struct id_list **un_current,fibreid *homepos,fibreid *currentpos,long *patrolregion,int nx,int ny)
{

  long flag,v1,v2,v3;  //variables to contain data from zone map
  int idnum; //id number from homepos
  int index,ind;  //index in homepos and currentpos respectively

  struct id_list *ptr_c,*ptr_f;

  int ncount,ii,jj,i;

    //set temporary pointer to the unassigned list

  ptr_f=(*un_fibre);

    //cycle through each of the home positions.   
    //if it's number is only assignable to one 
    //measured position, we know which it is. 

  while (ptr_f != NULL)
    {
      
      //get the ID number and index of the current COBRA

      idnum=ptr_f->idnum;
      index=ptr_f->index;

      ncount=0;

      //cycle through measured positions and look for possible matches for that fibre
      ptr_c=(*un_current);
      while(ptr_c != NULL)
	{
	  	  
	  //get position on map
	  ii=(int)(currentpos[ptr_c->index].x);
	  jj=(int)(currentpos[ptr_c->index].y);

	  //get the position 
	  flag=patrolregion[ii*ny+jj];

	  //parse the value to get the potential COBRA matches for the position
	  v1=flag % 10000;
	  v2 = ((flag - v1) % 100000000)/10000;
	  v3=flag/100000000;
	      

	  //is it a potential match for this 
	  if((idnum==v1)||(idnum==v2)||(idnum==v3))
	    {
	      //if the match is unambiguous (this is the only id that works for this
	      //measured position, break out of the loop
	      if(flag < 3000)
		{
		  ind=ptr_c->index;

		  ptr_c=(*un_current);
		  ncount=-1;
		  break;
		}
	      //otherwise, add the homepos id to the list of potential matches for
	      //the current position
	      else
		{

		  ncount+=1;
		  ind=ptr_c->index;
		  add_node_list(&ptr_c,(int)v1);
		  add_node_list(&ptr_c,(int)v2);
		  if(flag > 300000000)
		    {
		      add_node_list(&ptr_c,(int)v3);
		    }
		  add_node_list(&ptr_f,ind);
		 
		}

	    }
	  ptr_c=ptr_c->next;      
	}

      //incrment ptr_f first so we don't deallocate it

      if(ncount==-1)
	{

	  ptr_f=ptr_f->next;

		  
	  //assign the ID number
	  currentpos[ind].idnum=idnum;
	  update_list(as_fibre,un_fibre,un_current,idnum,index,ind,ind,homepos,currentpos);
	}
      else
	{
	  ptr_f=ptr_f->next;

	}

      //next fibre in list
    }

  //now we have identified those fibres which only have one possible
  //identification. 

  return;
}
/*--------------------------------------------------------*/

void  newly_singled(struct id_list **as_fibre,struct id_list **un_fibre,struct id_list **un_current,fibreid *homepos,fibreid *currentpos,int nx,int ny)
{

  int idnum; //id number from homepos
  int index,ind;  //index in homepos and currentpos respectively
  int ii;

  struct id_list *ptr_f,*ptr_t;
  struct match_list *ptr_a;

  int np,check;

  int newfit=1;
  int remove;
  //we need to loop this to account for newly singled ids

  while(newfit==1)
    {
      newfit=0;
      ptr_f=(*un_fibre);
      while (ptr_f != NULL)
	{
	  idnum=ptr_f->idnum;
	  index=ptr_f->index;

	  ptr_a=ptr_f->match;
	  np=0;
	  while(ptr_a != NULL)
	    {
	      np=np+1;
	      ptr_a=ptr_a->next;
	    }
	  check=0;
	  //has only one possible match
	  if(np==1)
	    {
	      
	      newfit=1;
	      index=ptr_f->index;
	      idnum=ptr_f->idnum;
	      ind=ptr_f->match->idnum;
	      
	      remove=0;
	      ptr_f=ptr_f->next;

	      update_list(as_fibre,un_fibre,un_current,idnum,index,ind,ind,homepos,currentpos);
	      
	      ptr_t=(*un_fibre);
	      while(ptr_t != NULL)
		{
		  ii=remove_node_list(&ptr_t->match,ind);
		  ptr_t=ptr_t->next;
		}
	      
	    }
	  else
	    {
	      ptr_f=ptr_f->next;
	    }
	}
    }

  return;

}
/*--------------------------------------------------------*/

void  nearest_point(struct id_list **as_fibre,struct id_list **un_fibre,struct id_list **un_current,fibreid *homepos,fibreid *currentpos,int nx,int ny)
{

  int idnum; //id number from homepos
  int index,ind;  //index in homepos and currentpos respectively
  int indm; //index of minimum distance

  struct id_list *ptr_f;
  struct match_list *ptr_a;

  double xx,yy,xc,yc;

  double dd,ddm; //distance calculations

  ptr_f=(*un_fibre);

  //now we're down to a list of fibres that have more than one match - 
  //for the fine move, we look for the closest previous position. 
  while (ptr_f != NULL)
    {
      idnum=ptr_f->idnum;
      index=ptr_f->index;

      //get the previous position for the fibre

      xc=homepos[index].xp;
      yc=homepos[index].yp;


      //now cycle through the potential matches
      ptr_a=ptr_f->match;
      ddm=500.;

      while (ptr_a != NULL)
	{
	  //current position of fibre
	  ind=ptr_a->idnum;
	  xx=currentpos[ind].x;
	  yy=currentpos[ind].y;
	  //find the nearest previous point
	  dd=sqrt((xx-xc)*(xx-xc)+(yy-yc)*(yy-yc));
	  //we want the closest position
	  if(dd<ddm)
	    {
	      ddm=dd;
	      indm=ind;
	    }

	  ptr_a=ptr_a->next;
	}

      ptr_f=ptr_f->next;

      update_list(as_fibre,un_fibre,un_current,idnum,index,ind,indm,homepos,currentpos);
      	      
    }
  //and repeat for each as yet unassigned fibre
}



int *matchrot(struct fibreid *homepos,struct fibreid *curpos,double theta,double xc,double yc,double xt,double yt,int nhome, int ncur)
{

  /* program to match home positions and identify fixed fiducial
     fibres in a field that is rotated around a point.

     homepos (fibreid): expected centre positions, without field rotation,

     curpos (fibreid): measured centre positions, with field rotation,

     xc,yc - centre of rotation, in pixels 

     theta - approximate angle of rotation about xc, yc
   */

  int i,j;  //loop counter

  float dd,dd1; //distances measurements (in pixels)

  //first, rotate the home positions to the current orientation

  //new positions
  double *xn=malloc(nhome*sizeof(double));
  double *yn=malloc(nhome*sizeof(double));

  //flags, for poor matches
  //int *ids=malloc(nhome*sizeof(int));
  int *dflag=malloc(nhome*sizeof(int));
  
  //flag for non-matching distance (in pixels)
  int df=2;

  //calculate the expected position  with field rotation
  for (i=0;i<nhome;i++)
    {
      xn[i]=(homepos[i].x-xc)*cos(theta)-(homepos[i].y-yc)*sin(theta)+xc+xt;
      yn[i]=(homepos[i].x-xc)*sin(theta)+(homepos[i].y-yc)*cos(theta)+yc+yt;
    }

  //Find the nearest point for each position
   for (i=0;i<nhome;i++)
     {
       dd=1000;
       for (j=0;j<ncur;j++)
	 {
	   dd1=(xn[i]-curpos[i].x)*(xn[i]-curpos[i].x)+(yn[i]-curpos[i].y)*(yn[i]-curpos[i].y);
	   if (dd1 < dd)
	     {
	       dd=dd1;
	       curpos[i].idnum=homepos[j].idnum;
	     }
	 }

       //flag anomalous nearest distances
       if (dd > df)
	 {
	   dflag[i]=dd;
	 }

     }
   return 0;
}

/*------------------------------------------------------------------------*/

int fibre_identify_fine(struct fibreid *homepos,struct fibreid *currentpos,long *patrolregion, int nx, int ny,int ncur,int nhome)
{

  /*

    routine to identify fibres for later (smaller) moves.

    INPUT

    homepos - home positions of the fibres, with id numbers, rotated as needed
    currentpos - measured fibre id positions
    patrolregion - a map of potential areas where confusion is possible

    OUTPUT
    
    fibreid *newid - identified positions

  */

  int i; // loop variables, pixel values


  //pointers to lists of fibre or centroid ids

  struct id_list *un_current=NULL;
  struct id_list *as_fibre=NULL;
  struct id_list *un_fibre=NULL;

  //create a linked list of the ID numbers and their respective
  //indices in the homepos array. this is the unassigned fibre list.

  for (i=0;i<nhome;i++)
    {
      add_node(&un_fibre,homepos[i].idnum,i);
    }

  //create a corresponding list of the unassigned measured positions

  for (i=0;i<ncur;i++)
    {
      add_node_cur(&un_current,i);      
    }


  //For fine motion, we assume that fibres have been moved out of 
  //hidden positions. 

  initial_position(&as_fibre,&un_fibre,&un_current,homepos,currentpos,patrolregion,nx,ny);

  newly_singled(&as_fibre,&un_fibre,&un_current,homepos,currentpos,nx,ny);
  nearest_point(&as_fibre,&un_fibre,&un_current,homepos,currentpos,nx,ny);

  return 0;
}


/*--------------------------------------------------------*/

int fibre_identify_coarse(struct fibreid *homepos,struct fibreid *currentpos,long *patrolregion,int *arc_image,int *spiral,int *adjacent_fibre_list,int npval,int nx,int ny,int ncur,int nhome)
{

  /*

    routine to identify fibres for the first (and largest) move. 

    INPUT

    homepos - home positions of the fibres, with id numbers, rotated as needed
    currentpos - measured fibre id positions
    patrolregion - a map of potential areas where confusion is possible

    OUTPUT
    
    fibreid newid - identified positions

  */

  //pointers to lists of fibre or centroid ids

  int i;
  struct id_list *un_current=NULL;
  struct id_list *as_fibre=NULL;
  struct id_list *un_fibre=NULL;
  

  //create a linked list of the ID numbers and their respective
  //indices in the homepos array. this is the unassigned fibre list.

  //printf("Entering identify_fibre_coarse\n");
  //printf("Setting Up Lists\n");

  for (i=0;i<nhome;i++)
    {
      add_node(&un_fibre,homepos[i].idnum,i);
    }

  //create a corresponding list of the unassigned measured positions

  for (i=0;i<ncur;i++)
    {
      add_node_cur(&un_current,i);      
    }
  


  initial_position(&as_fibre,&un_fibre,&un_current,homepos,currentpos,patrolregion,nx,ny);

  //i=print_list(un_fibre);
  zero_matches(&as_fibre,&un_fibre,&un_current,homepos,currentpos,arc_image,nx,ny);

  //i=print_list(&un_fibre);

  struct id_list *ptr_f;
  ptr_f=un_fibre;

  newly_singled_coarse(&as_fibre,&un_fibre,&un_current,homepos,currentpos,arc_image,nx,ny,adjacent_fibre_list);
  confused_fibres(&as_fibre,&un_fibre,&un_current,homepos,currentpos,arc_image,nx,ny,spiral,npval);

  return 0;
}

/*--------------------------------------------------------------------*/

void  zero_matches(struct id_list **as_fibre,struct id_list **un_fibre,struct id_list **un_current,fibreid *homepos,fibreid *currentpos,int *arc_image,int nx,int ny)

{  

  int idnum; //id number from homepos
  int index;  //index in homepos and currentpos respectively

  struct id_list *ptr_f;

  double xc,yc;

  int x_mask_offset=20;
  int y_mask_offset=20;

  ptr_f=(*un_fibre);

  while (ptr_f != NULL)
    {
      idnum=ptr_f->idnum;
      index=ptr_f->index;

      //no potential matches
      if(ptr_f->match==NULL)
	{	    

	  //printf("vc  %lf %lf %d %d\n",homepos[index].x,homepos[index].y,index,idnum);
	  ptr_f=ptr_f->next;

	  xc=homepos[index].x+x_mask_offset;
	  yc=homepos[index].y+y_mask_offset;

	  update_list_hidden(as_fibre,un_fibre,un_current,idnum,index,xc,yc,homepos,currentpos);

	  //check=is_hidden(arc_image,nx,ny,iyc,ixc,6,11,300,5);
	  //check=0;

	}
      else
	{
	  ptr_f=ptr_f->next;
	}
    }
}
/*--------------------------------------------------------*/

void  newly_singled_coarse_withmask(struct id_list **as_fibre,struct id_list **un_fibre,struct id_list **un_current,fibreid *homepos,fibreid *currentpos,int *arc_image,int nx,int ny,int *adjacent_fibre_list)


{

  int idnum; //id number from homepos
  int index,ind;  //index in homepos and currentpos respectively

  int check;
  struct id_list *ptr_f,*ptr_t;
  struct match_list *ptr_a;

  int ii;
  double xc,yc;
  int check1,remove,np;

  int x_mask_offset=20;
  int y_mask_offset=20;

  int newfit=1;
  //we need to loop this to account for newly singled ids

  while(newfit==1)
    {

      //count the number of possible matches
      newfit=0;
      ptr_f=(*un_fibre);
      while (ptr_f != NULL)
	{
	  idnum=ptr_f->idnum;
	  index=ptr_f->index;
	  
	  ptr_a=ptr_f->match;
	  np=0;
	  while(ptr_a != NULL)
	    {
	      np=np+1;
	      ptr_a=ptr_a->next;
	    }
	  check=0;
	  //has only one possible match
	  if(np==1)
	    {
	      
	      newfit=1;
	      index=ptr_f->index;
	      idnum=ptr_f->idnum;
	      ind=ptr_f->match->idnum;
	      
	      remove=0;
	      //check to see if adjacent fibres are assigned
	      check=adjacent_assigned((*un_fibre),index,adjacent_fibre_list);
	      //if the adjacent fibres are all assigned, then make the assignment
	      if(check==1)
		{

		  ptr_f=ptr_f->next;

		  update_list(as_fibre,un_fibre,un_current,idnum,index,ind,ind,homepos,currentpos);

		  remove=1;
		  
		}
	      //if they're not, check to see if the fibre hidden
	      else 
		{
		  xc=homepos[index].x+x_mask_offset;
		  yc=homepos[index].y+y_mask_offset;
		  check1=is_hidden(arc_image,nx,ny,(int)xc,(int)yc,8,11,250,5);
		  check1=0;
		  
	      //if the fibre is hidden, assign the mask position
		  if(check1==1)
		    {

		      ptr_f=ptr_f->next;

		      update_list_hidden(as_fibre,un_fibre,un_current,idnum,index,xc,yc,homepos,currentpos);
		      
		    }
		  //if it's not hidden, then we can assign the position
		  else
		    {

		      ptr_f=ptr_f->next;

		      update_list(as_fibre,un_fibre,un_current,idnum,index,ind,ind,homepos,currentpos);

		      remove=1;
		      
		    }
		}
	      
	      //IF A FIBRE HAS BEEN ASSIGNED, remove it from the possible fibres for
	      //the other unassigned fibres
	      if(remove==1)
		{
		  ptr_t=(*un_fibre);
		  while(ptr_t != NULL)
		    {
		      ii=remove_node_list(&ptr_t->match,ind);
		      ptr_t=ptr_t->next;
		    }
		  
		  
		}
	      
	    }
	  else
	    {
	      ptr_f=ptr_f->next;
	    }
	}
    }

}
/*--------------------------------------------------------------------*/

void  newly_singled_coarse(struct id_list **as_fibre,struct id_list **un_fibre,struct id_list **un_current,fibreid *homepos,fibreid *currentpos,int *arc_image,int nx,int ny,int *adjacent_fibre_list)
{

  int idnum; //id number from homepos
  int index,ind;  //index in homepos and currentpos respectively

  int check;
  struct id_list *ptr_f,*ptr_t,*ptr_g;
  struct match_list *ptr_a;

  int ii;
  double xc,yc;
  int check1,remove,np;

  int x_mask_offset=20;
  int y_mask_offset=20;

  int newfit=1;
  //we need to loop this to account for newly singled ids

  while(newfit==1)
    {

      //count the number of possible matches
      newfit=0;
      ptr_f=(*un_fibre);
      while (ptr_f != NULL)
	{

	  
	  idnum=ptr_f->idnum;
	  index=ptr_f->index;
	  ptr_a=ptr_f->match;

	  np=0;
	  while(ptr_a != NULL)
	    {
	      np=np+1;
	      ptr_a=ptr_a->next;
	      
	    }

	  check=0;
	  //has only one possible match
	  if(np==1)
	    {
	      
	      newfit=1;
	      index=ptr_f->index;
	      idnum=ptr_f->idnum;
	      ind=ptr_f->match->idnum;
	      
	      ptr_f=ptr_f->next;
	      update_list(as_fibre,un_fibre,un_current,idnum,index,ind,ind,homepos,currentpos);
	      remove=1;
	    	      
	      //IF A FIBRE HAS BEEN ASSIGNED, remove it from the possible fibres for
	      //the other unassigned fibres
	      if(remove==1)
		{

		  ptr_t=(*un_fibre);
		  while(ptr_t != NULL)
		    {
		      ii=remove_node_list(&ptr_t->match,ind);
		      ptr_t=ptr_t->next;
		    }
				  
		  
		}
	      
	    }
	  else
	    {
	      ptr_f=ptr_f->next;
	    }

	}
    }
}
/*--------------------------------------------------------------------*/

void  confused_fibres(struct id_list **as_fibre,struct id_list **un_fibre,struct id_list **un_current,fibreid *homepos,fibreid *currentpos,int *arc_image,int nx,int ny,int *spiral,int npval)


{

  int idnum; //id number from homepos
  int index,ind;  //index in homepos and currentpos respectively
  int indm;

  int *pos;
  struct id_list *ptr_f;
  struct match_list *ptr_a;

  double xc,yc,xx,yy;

  double dd,ddm; //distance calculations

  //now we're down to a list of fibres that have more than one match - 
  //we need to check where that fibre came from. 
  ptr_f=(*un_fibre);
  while (ptr_f != NULL)
    {

      //get the centre position for the fibre
      index=ptr_f->index;
      idnum=ptr_f->idnum;
      xc=homepos[index].x;
      yc=homepos[index].y;
      
      //now cycle through the potential matches
      ptr_a=ptr_f->match;
      ddm=500.;

      while (ptr_a != NULL)
	{
	  //current position of fibre
	  
	  ind=ptr_a->idnum;
	  xx=currentpos[ind].x;
	  yy=currentpos[ind].y;

	  //find the direction from which it came
	  pos=wherefrom(arc_image,nx,ny,yy,xx,spiral,npval);
	  //printf("  fibre position=%ld %lf %lf\n",ind,xx,yy);
	  //printf("   wherefrom=%d %d %d\n",ind,pos[0],pos[1]);
	  
	  //distance betwen the direction and the centre position
	  dd=sqrt((xc-pos[1])*(xc-pos[1])+(yc-pos[0])*(yc-pos[0]));
	  
	  //printf("   dd=%lf\n",dd);

	  //we want the closest position
	  if(dd<ddm)
	    {
	      ddm=dd;
	      indm=ind;
	    }


	  ptr_a=ptr_a->next;
	}

      //assignment and list management

      ptr_f=ptr_f->next;
      update_list(as_fibre,un_fibre,un_current,idnum,index,indm,indm,homepos,currentpos);

    }
}
