#include <math.h>
#include <iostream>
#include "nr.h"
#include <stdlib.h>
#include "../solveqp.h"

using namespace std;

long geo(long in)
//geo: greater/equal one
{
  if(in>1)
    return(in);
  else
    return(1);
}

long len(long in,long n)
//gen: less/equal n
{
  if(in<n)
    return(in);
  else
    return(n);
}

void choldc(double **a, int n)
{
  int i,j,k;
  double sum=0;
  double *p = dvector(1,n);

  for (i=1;i<=n;i++) 
    {
    for (j=i;j<=n;j++) 
      {
	for (sum=a[i][j],k=i-1;k>=1;k--)
	  sum -= a[i][k]*a[j][k];
	if (i == j) {
	  if (sum <= 0.0)
	    cerr << "choldc failed --> matrix not positive definite\n";
	  p[i]=sqrt(sum);
	} 
	else
	  {
	    a[j][i]=sum/p[i];
	  }
      }
    }


  for (i=1;i<=n;i++)
    a[i][i] = p[i];

  for (i=1;i<=n;i++) 
    {
      for (j=i+1;j<=n;j++)
	a[i][j]=0;
    }
  free_dvector(p,1,n);
}

//Modifizierte Cholesky-Zerlegung fuer Bandmatrizen (BW=5)

void Choldc(Bandmatrix *a, int n, double p[])
{
  long k,l,m;
  double sum;

  for(k=1;k<=n;k++)
    {
      for(l=k;l<=len(k+2,n);l++)
	{
	  for(sum=a->r(k,l),m=k-1;m>=1;m--)
	    {
	      if(abs(k-m)<=2 && abs(l-m)<=2)
		sum-=(a->r(k,m)*a->r(l,m));
	    }
	  if(k==l)
	    {
	      if(sum<=0.0)
		{
		  cerr << "Sum : " << sum << endl;
		  cerr << "choldc failed\n"; exit(1);
		}
	      p[k]=sqrt(sum);
	    }
	  else
	    a->w(sum/p[k],l,k);
	}
    }
}

