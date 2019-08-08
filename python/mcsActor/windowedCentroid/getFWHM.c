

int getFWHM(struct cand_point *pointList,double *fx,double *fy)
{

  struct cand_point *curr=pointList;

  int sumX=0,sumY=0,sumT=0,np=0;
  
  while (curr != NULL)
    {
      sumX+=curr->x2;
      sumY+=curr->y2;
      sumT+=curr->xy;
      np+=1;
      
    }

  double x2=(double)sumX/(double)np;
  double y2=(double)sumY/(double)np;
  double theta=(double)sumT/(double)np;

  (*fx)=x2*abs(sin(theta));
  (*fy)=y2*abs(cos(theta));

  return 0;
  
}
