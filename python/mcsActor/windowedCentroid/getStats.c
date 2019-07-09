


int getStats(int *image,int nx, int ny, double *mean, double *std)
{

  int i,j,pixSum=0,pixSum2=0;

  for (i=0;i<nx;i++)
    {
      for (j=0;j<ny;j++)
	{
	  pixInd=getInd2D(i,j,nx);
	  pixSum+=image[pixInd];
	  pixSum2+=image[pixInd]*image[pixInd];
	}
    }
  double np=(double)nx*(double)ny;
  (*mean)=(double)pixSum/np;
  (*std)=sqrt((pixSum2/np-(*mean)*(*mean))/np);
	  
  return 0;
  
}
