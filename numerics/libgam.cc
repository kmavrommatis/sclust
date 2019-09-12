#include<fstream>
#include<iostream>
#include<stdlib.h>
#include<math.h>
#include "nr.h"

using namespace std;

//----------------------------------------------------------------------
//                  SPEZIELLE KLASSE FUER BANDMATRIZEN
//----------------------------------------------------------------------


//Konstruktor
Bandmatrix::Bandmatrix(long ndata,long bandbreite) : n(ndata), bw(bandbreite)
{
  long k,l;
  if(bw<0)
    {
      cerr << "Bandwidth insufficient !!\n";
      exit(1);
    }
  band=dmatrix(1,n,0,2*bw+1);
  //die Bandmatrix mit Null initialisieren
  for(k=1;k<=n;k++)
    {
      for(l=0;l<=2*bw+1;l++)
	{
	  band[k][l]=0.;
	}
    }
}

//Lesen der Eintraege
double Bandmatrix::r(long i,long j)
{
  if(abs(i-j)>bw)
    {
//      cerr << "*";
      return(0.);
    }
  else
    {
      return(band[i][bw-(i-j)]);
    }
}

//Schreiben der Eintraege
void Bandmatrix::w(double value,long i,long j)
{
  if(abs(i-j)>bw);
  else
    {
      band[i][bw-(i-j)]=value;
    } 
}

//Destruktor
Bandmatrix::~Bandmatrix()
{
  free_dmatrix(band,1,n,0,bw);
}


//------------------------------------------------------------

void splines(double *x,double *y,double *w,long n,double alpha,double *g,double *gam)
//Berechnet smoothing splines (natural cubic splines) nach dem
//Reinsch-Algorithmus
// double *x : n-Vektor, design points
// double *y : n-Vektor, Daten
// double *w : n-Vektor, Gewichte
// long n : #Datenpunkte
// double alpha: smoothing-Parameter
// double *g : n-Vektor, Ausgabewerte
// double *gam: n-2 Vector, Ausgabe zweite Ableitung
{
  double *h,*b,*p,dummy;
  long k,l,m;

  if(n<=2)
    {
      cerr << "Error : n < 3 !!! \n";
      exit(1);
    }
  if(alpha<0)
    alpha=alpha*(-1);

  h=dvector(1,n-1);
  Bandmatrix Q(n,2),A(n-2,2);
   
  b=dvector(1,n-2);
  p=dvector(1,n-2);

  //Initialisierung von h,b
  //h
  for(k=1;k<=n-1;k++)
    h[k]=x[k+1]-x[k];
  //b
  for(k=1;k<=n-2;k++)
    b[k]=(y[k+2]-y[k+1])/h[k+1]-(y[k+1]-y[k])/h[k];
  // Q & A

  for(k=1;k<=n-2;k++)
    {
      Q.w(1./h[k],k,k);
      Q.w(-1./h[k]-1./h[k+1],k+1,k);
      Q.w(1./h[k+1],k+2,k);
    }

  // Ausnutzen der Bandstruktur fuer schnelle Berechnung von Q^t*Q

  for(l=-2;l<=2;l++)
    {
      for(m=-2;m<=2;m++)
	{
	  for(k=1;k<=n;k++)
	    {
	      if(k+l >=1 && k+m >=1 && k+l <=n-2 && k+l <=n-2)
		A.w(A.r(k+l,k+m)+Q.r(k,k+l)*(1./w[k])*Q.r(k,k+m),k+l,k+m);
	    }
	}
    }

  for(k=1;k<=n-2;k++)
    {
      for(l=-2;l<=2;l++)
	{
	  if(k+l >=1 && k+l <=n-2)
	  A.w(alpha*A.r(k,k+l),k,k+l);
	}
    }

  for(k=1;k<=n-2;k++)
    {
      A.w(A.r(k,k)+(1./3.)*(h[k]+h[k+1]),k,k);
      if(k<=n-3)
	{
	  A.w(A.r(k,k+1)+h[k+1]/6.,k,k+1);
	  A.w(A.r(k+1,k)+h[k+1]/6.,k+1,k);
	}
    }
  //Loesen der Gleichungen mittels Cholesky-Zerlegung
  Choldc(&A,n-2,p);
  Cholsl(&A,n-2,p,b,gam);
   
   //Berechen der Ausgabewerte g

  for(k=1;k<=n;k++)
    {
      dummy=0.;
      for(l=1;l<=n-2;l++)
	{
	  dummy+=Q.r(k,l)*gam[l];
	}
      g[k]=y[k]-alpha*(1./w[k])*dummy;
    }
  
   
   
   free_dvector(h,1,n-1);
   free_dvector(b,1,n-2);
   free_dvector(p,1,n-2);
}

long make_w(double *xs,double*ys,long *ar,long n,double *w)
{
  long k,l;
  double x;

  l=1;
  x=xs[1];
  ar[1]=1;
  for(k=1;k<=n;k++)
    w[k]=1.;
  for(k=2;k<=n;k++)
    {
      if(xs[k]==x)
	{
	  w[l]++;
	  ys[l]+=ys[k];
	}
      else
	{
	  l++;
	  x=xs[k];
	  ys[l]=ys[k];
	  xs[l]=xs[k];
	}
      ar[k]=l;
    }
  for(k=1;k<=l;k++)
    ys[k]*=1./(double)w[k];
  return(l);
}

double back_fit(double **x,double *y,long n,long order,double *alpha,double **g,double **gam,double **x_out,double **y_out,double **w,double **g_out,long *n_out,int flag)
  //Eine Back-Fitting Iteration, liefert die Summe ueber Residuenquadrat zurueck -> Abbruchbed.
  //GCV IST HIER IMPLEMENTIERT, Art des Backfitting : Gauss-Seidel
  // und benutzt natural qubic splines als Glaetter
  // double **x : order,n-Matrix, design points
  // double *y : n-Vektor, Daten ; n : #Datenpunkte
  // int order : Modell-Ordnung
  // double *alpha: order-Vektor smoothing-Parameter Vektor fuer jede Funktion
  // double **g : order,n-Matrix der Ausgabewerte - nur fuer init und innere Verarb.
  // doube **gam : order,n-2 - Matrix zweite Ableitung
  // double **x_out  : order,n-Matrix, design points der Ausgabe (sortiert)
  // double **y_out : order,n-Matrix liefert die "Datenverteilung" wieder
  // double **w : order,n-Matrix die Gewichte fuer die Glaettung
  // double **g_out : order,n-Matrix, Ausgabewerte (sortiert)
  // long **n_out : order-Vektor #Datenpunkte der versch. Funktionen
  // int flag : Ist flag = 1 -> Backfit. ohne GCV ; ist flag = 2 -> Backfit. mit GCV
{
  long k,l,m,*indx,*ar;
  double *yj,*ys,*gj,*gamj,*xj,dummy,res;
  double *xs,*ws,alpha_zw;
  yj=dvector(1,n);
  ys=dvector(1,n);
  gj=dvector(1,n);
  gamj=dvector(1,n-2);
  xj=dvector(1,n);
  xs=dvector(1,n);
  ws=dvector(1,n);
  indx=lvector(1,n);
  ar=lvector(1,n);

  for(k=1;k<=order;k++)
    {
      for(m=1;m<=n;m++)
	yj[m]=y[m];
      for(l=1;l<=order;l++)
	{
	  if(l!=k)
	    {
	      for(m=1;m<=n;m++)
		yj[m]-=g[l][m];
	    }
	  else
	    {
	      for(m=1;m<=n;m++)
		xj[m]=x[k][m];
	    }
	}
      indexx(n,xj,indx);  // sorting
      for(l=1;l<=n;l++)
	{
	  ys[l]=yj[indx[l]];
	  xs[l]=xj[indx[l]];
	}
      n_out[k]=make_w(xs,ys,ar,n,ws); //Gewichtung : Tied design points !!
      if(flag==1)
	splines(xs,ys,ws,n_out[k],alpha[k],gj,gamj); 
      else if(flag==2)
	{
	  alpha_zw=alpha[k];
	  splines_gcv(xs,ys,ws,n_out[k],&alpha_zw,gj,gamj);
	  alpha[k]=alpha_zw;
	}
      else
	{
	  cerr << "Error in flag !!!\n";
	  exit(1);
	}

      for(m=1;m<=n_out[k];m++)
	{
	  x_out[k][m]=xs[m];
	  w[k][m]=ws[m];
	  y_out[k][m]=ys[m];
	  g_out[k][m]=gj[m];
	  if(m<=n_out[k]-2)
	    gam[k][m]=gamj[m];
	}
      for(m=1;m<=n;m++)
	g[k][indx[m]]=gj[ar[m]];
    }

  res=0.;
  for(k=1;k<=n;k++)
    {
      dummy=0.;
      for(l=1;l<=order;l++)
	{
	  dummy+=g[l][k];
	}
      res+=pow(y[k]-dummy,2);
    }
  free_dvector(yj,1,n);
  free_dvector(ys,1,n);
  free_dvector(ws,1,n);
  free_dvector(gj,1,n);
  free_dvector(gamj,1,n-2);
  free_dvector(xj,1,n);
  free_dvector(xs,1,n);
  free_lvector(indx,1,n);

  return(res);
}

double plot_splines(double *x,double *g,double *gam,long n,double t)
// Auswerten der Splines
// double *x : n-Vektor Design-Points
// double *g : n-Vektor Stuetzstellen
// double *gam : n-Vektor zweite Ableitung Vorsicht: Andere Notation gam[1]=0 und gam[n]=0
// long n : #Punkte
// double t : Auswertstelle
{
  double out,dum;
  long k;

  if(t<=x[1])
    {
      dum=(g[2]-g[1])/(x[2]-x[1])-(1./6.)*gam[2]*(x[2]-x[1]);
      out=g[1]-(x[1]-t)*dum;
    }
  else if(t>=x[n])
    {
      dum=(g[n]-g[n-1])/(x[n]-x[n-1])+(1./6.)*gam[n-1]*(x[n]-x[n-1]);
      out=g[n]+(t-x[n])*dum;
    }
  else
    {
      k=1;
      while(t>x[k])
	{
	  k++;
	}
      k--;
      dum=x[k+1]-x[k];
      out=((t-x[k])*g[k+1]+(x[k+1]-t)*g[k])/dum-(1./6.)*(t-x[k])*(x[k+1]-t)*((1+(t-x[k])/dum)*gam[k+1]+(1+(x[k+1]-t)/dum)*gam[k]);
    }
  return(out);
}

double deriv_splines(double *x,double *g,double *gam,long n,double t)
// Auswerten der Ableitung der Splines
// double *x : n-Vektor Design-Points
// double *g : n-Vektor Stuetzstellen
// double *gam : n-Vektor zweite Ableitung Vorsicht: Andere Notation gam[1]=0 und gam[n]=0
// long n : #Punkte
// double t : Auswertstelle
{
  double out,dum,a,b,c,d;
  long k;

  if(t<=x[1])
    {
      dum=(g[2]-g[1])/(x[2]-x[1])-(1./6.)*gam[2]*(x[2]-x[1]);
      out=-dum;
    }
  else if(t>=x[n])
    {
      dum=(g[n]-g[n-1])/(x[n]-x[n-1])+(1./6.)*gam[n-1]*(x[n]-x[n-1]);
      out=dum;
    }
  else
    {
      k=1;
      while(t>x[k])
        {
          k++;
        }
      k--;
      dum=x[k+1]-x[k];
      a=(t-x[k]);
      b=(x[k+1]-t);
      c=((1+(t-x[k])/dum)*gam[k+1]+(1+(x[k+1]-t)/dum)*gam[k]);
      d=gam[k+1]/dum-gam[k]/dum;
      out=(g[k+1]-g[k])/dum-(1./6.)*(b*c-a*c+a*b*d);
    }
  return(out);
}


// Nochmals Backfitting aber diesmal mit Kernschaetzer

double K(double u,double h)
{
  double out;
  if(fabs(u)<=h)
    {
      if(u<=0)
	out=1/h*(1./h*u+1);
      else
	out=1/h*(-1./h*u+1);
    }
  else
    out=0.;
  return(out);
}

void kernel_smooth(double *x,double *y,long n,double h,double *w,double *g)
{
  double nenner,zaehler;
  long k,l,bound;
  
  for(k=1;k<=n;k++)
    {
      nenner=K(x[k]-x[k],h)*w[k];
      zaehler=K(x[k]-x[k],h)*w[k]*y[k];

      for(l=1;l<=n;l++)
	{
	  if(k+l>n)
	    break;
	  if((x[k+l]-x[k])> h)
	    break;
	  nenner+=K(x[k]-x[k+l],h)*w[k+l];
	  zaehler+=K(x[k]-x[k+l],h)*w[k+l]*y[k+l];
	}
      for(l=1;l<=n;l++)
	{	  
	  if(k-l<=0)
	    break;
	  if((x[k]-x[k-l])>h )
	    break;
	  nenner+=K(x[k]-x[k-l],h)*w[k-l];
	  zaehler+=K(x[k]-x[k-l],h)*w[k-l]*y[k-l];
	}
      g[k]=zaehler/nenner;
    }
  
} 

double back_fit_kernel(double **x,double *y,long n,long order,double *h,double **g,double **x_out,double **y_out,double **w,double **g_out,long *n_out)
  //Eine Back-Fitting Iteration, liefert die Summe ueber Residuenquadrat zurueck -> Abbruchbed.
  // und benutzt Kernschaetzer als Glaetter
  // double **x : order,n-Matrix, design points
  // double *y : n-Vektor, Daten ; n : #Datenpunkte
  // int order : Modell-Ordnung
  // double *h: order-Vektor smoothing-Parameter Vektor fuer jede Funktion
  // double **g : order,n-Matrix der Ausgabewerte - nur fuer init und innere Verarb.
  // double **x_out  : order,n-Matrix, design points der Ausgabe (sortiert)
  // double **y_out : order,n-Matrix liefert die "Datenverteilung" wieder
  // double **w : order,n-Matrix die Gewichte fuer die Glaettung
  // double **g_out : order,n-Matrix, Ausgabewerte (sortiert)
  // long **n_out : order-Vektor #Datenpunkte der versch. Funktionen
{
  long k,l,m,*indx,*ar;
  double *yj,*ys,*gj,*xj,dummy,res;
  double *xs,*ws;

  yj=dvector(1,n);
  ys=dvector(1,n);
  gj=dvector(1,n);
  xj=dvector(1,n);
  xs=dvector(1,n);
  ws=dvector(1,n);
  indx=lvector(1,n);
  ar=lvector(1,n);

  for(k=1;k<=order;k++)
    {
      for(m=1;m<=n;m++)
	yj[m]=y[m];
      for(l=1;l<=order;l++)
	{
	  if(l!=k)
	    {
	      for(m=1;m<=n;m++)
		yj[m]-=g[l][m];
	    }
	  else
	    {
	      for(m=1;m<=n;m++)
		xj[m]=x[k][m];
	    }
	}
      indexx(n,xj,indx);  // sorting
      for(l=1;l<=n;l++)
	{
	  ys[l]=yj[indx[l]];
	  xs[l]=xj[indx[l]];
	}
      n_out[k]=make_w(xs,ys,ar,n,ws);
      kernel_smooth(xs,ys,n_out[k],h[k],ws,gj); 
      for(m=1;m<=n_out[k];m++)
	{
	  x_out[k][m]=xs[m];
	  w[k][m]=ws[m];
	  y_out[k][m]=ys[m];
	  g_out[k][m]=gj[m];
	}
      for(m=1;m<=n;m++)
	g[k][indx[m]]=gj[ar[m]];
    }

  res=0.;
  for(k=1;k<=n;k++)
    {
      dummy=0.;
      for(l=1;l<=order;l++)
	{
	  dummy+=g[l][k];
	}
      res+=pow(y[k]-dummy,2);
    }

  free_dvector(yj,1,n);
  free_dvector(ys,1,n);
  free_dvector(ws,1,n);
  free_dvector(gj,1,n);
  free_dvector(xj,1,n);
  free_dvector(xs,1,n);

  return(res);
}


double Score(double alpha,double *x,double *y,double *w,double *b,double *g,double *gam,Bandmatrix *Q,Bandmatrix *R,long n)
  //Berechnet wird der GCV-Score
{
  double score=0.,dummy;
  double spur=0.; // tr(S)
  long k,l,m;
  double *d;

  Bandmatrix A(n-2,2),B(n-2,2);
  d=dvector(1,n-2);

  //Spline berechnen (Reinsch Algorithmus)

  //Uebertragen der Elemente von R auf A
  for(k=1;k<=n-2;k++)
    {
      for(l=-2;l<=2;l++)
	{
	  if(k+l>=1 && k+l <=n-2)
	    A.w(R->r(k,k+l),k,k+l);
	}
    }
  
  for(l=-2;l<=2;l++)
    {
      for(m=-2;m<=2;m++)
	{
	  for(k=1;k<=n;k++)
	    {
	      if(k+l >=1 && k+m >=1 && k+l <=n-2 && k+l <=n-2)
		A.w(A.r(k+l,k+m)+alpha*Q->r(k,k+l)*(1./w[k])*Q->r(k,k+m),k+l,k+m);
	    }
	}
    }

  Choldc(&A,n-2,d);
  Cholsl(&A,n-2,d,b,gam);  

  
  for(k=1;k<=n;k++)
    {
      g[k]=y[k];
      for(l=k-2;l<=k+2;l++)
	{
	  if(l>=1 && l<=n-2)
	    g[k]-=alpha*(1./w[k])*Q->r(k,l)*gam[l];
	}
    }
  //Ende Spline berechnen
  //cerr << "*";
  //Spur berechnen nach nach Hutchison & de Hoog
  B.w(1./(d[n-2]*d[n-2]),n-2,n-2);
  B.w((-A.r(n-2,n-3)*B.r(n-2,n-2))/d[n-3],n-3,n-2);
  B.w(B.r(n-3,n-2),n-2,n-3);
  B.w(((1./d[n-3])-A.r(n-2,n-3)*B.r(n-3,n-2))/d[n-3],n-3,n-3);

  for(k=n-4;k>=1;k--)
    {
      B.w((-A.r(k+1,k)*B.r(k+1,k+2)-A.r(k+2,k)*B.r(k+2,k+2))/d[k],k,k+2);
      B.w(B.r(k,k+2),k+2,k);
      B.w((-A.r(k+1,k)*B.r(k+1,k+1)-A.r(k+2,k)*B.r(k+1,k+2))/d[k],k,k+1);
      B.w(B.r(k,k+1),k+1,k);
      B.w(((1./d[k])-A.r(k+1,k)*B.r(k,k+1)-A.r(k+2,k)*B.r(k,k+2))/d[k],k,k);
    }


  
  for(k=1;k<=n;k++)
    {
      dummy=0.;
      for(l=-2;l<=2;l++)
	{
	  for(m=-2;m<=2;m++)
	    {
	      if(k+l >=1 && k+m >=1 && k+l <=n-2 && k+l <=n-2)
		dummy+=Q->r(k,k+l)*B.r(k+l,k+m)*Q->r(k,k+m);
	    }
	}
      spur+=1.+alpha*(1./w[k])*dummy;
    }

  free_dvector(d,1,n-2);
 
  //Ende Spur berechnen

  //Score berechnen
  for(k=1;k<=n;k++)
    {
      score+=w[k]*pow(y[k]-g[k],2);
    }
  score*=1./((double)n*pow(1.-spur/(double)n,2));
  return(score);
  
}

//Golden Section Minimierungsroutiene -> Recipes

#define RR 0.61803399
#define C (1.0-RR)
#define SHFT2(a,b,c) (a)=(b);(b)=(c);
#define SHFT3(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

double golden(double ax, double bx, double cx, double tol,double *x,double *y,double *w,double *b,double *g,double *gam,Bandmatrix *Q,Bandmatrix *R,long n)
{
  double f1,f2,x0,x1,x2,x3;
  
  x0=ax;
  x3=cx;
  if (fabs(cx-bx) > fabs(bx-ax)) {
    x1=bx;
    x2=bx+C*(cx-bx);
  } else {
    x2=bx;
    x1=bx-C*(bx-ax);
  }
  f1=log(Score(x1,x,y,w,b,g,gam,Q,R,n));
  f2=log(Score(x2,x,y,w,b,g,gam,Q,R,n));
  while (fabs(x3-x0) > tol*(fabs(x1)+fabs(x2))) {
    if (f2 < f1) {
      SHFT3(x0,x1,x2,RR*x1+C*x3)
	SHFT2(f1,f2,log(Score(x2,x,y,w,b,g,gam,Q,R,n)))
	} else {
	  SHFT3(x3,x2,x1,RR*x2+C*x0)
	    SHFT2(f2,f1,log(Score(x1,x,y,w,b,g,gam,Q,R,n)))
	    }
  }
  if (f1 < f2) {
    
    return x1;
  } else {
    
    return x2;
  }
}

#undef C
#undef RR
#undef SHFT2
#undef SHFT3


void splines_gcv(double *x,double *y,double *w,long n,double *alpha,double *g,double *gam)
//Berechnet smoothing splines (natural cubic splines) nach dem
//Reinsch-Algorithmus der Glaettungsparameter wird durch Generalized Crossvalidation
//geschaetzt.
// double *x : n-Vektor, design points
// double *y : n-Vektor, Daten
// double *w : n-Vektor, Gewichte
// long n : #Datenpunkte
// double alpha: smoothing-Parameter
// double *g : n-Vektor, Ausgabewerte
// double *gam: n-2 Vector, Ausgabe zweite Ableitung
{
  double *h,*b,dummy;


  long k,l,m;

  if(n<=2)
    {
      cerr << "Error : n < 3 !!! \n";
      exit(1);
    }
  h=dvector(1,n-1);
  b=dvector(1,n-2);
  Bandmatrix Q(n,2),R(n,2);
  //Initialisierung von h,b
  //h
  for(k=1;k<=n-1;k++)
    h[k]=x[k+1]-x[k];
  //b
  for(k=1;k<=n-2;k++)
    b[k]=(y[k+2]-y[k+1])/h[k+1]-(y[k+1]-y[k])/h[k];

  //Q
  for(k=1;k<=n-2;k++)
    {
      Q.w(1./h[k],k,k);
      Q.w(-1./h[k]-1./h[k+1],k+1,k);
      Q.w(1./h[k+1],k+2,k);
    }
  //R
  for(k=1;k<=n-2;k++)
    {
      R.w((1./3.)*(h[k]+h[k+1]),k,k);
      if(k<=n-3)
	{
	  R.w((1./6.)*h[k+1],k,k+1);
	  R.w((1./6.)*h[k+1],k+1,k);
	}
    }
  
  
  //Minimum vom Score berechnen
  double ax=1e-10;
  double cx=1e10;

  //cerr << "GCV : ";
  
  (*alpha)=golden(ax,(*alpha),cx,1e-4,x,y,w,b,g,gam,&Q,&R,n);
  
  //Golden Section Minimierungsroutine

  //cerr << "\n Alpha = " << *alpha << endl;
  free_dvector(h,1,n-1);
  free_dvector(b,1,n-2);
}



void error_band(double **x,long n,long order,double *alpha,double **g,double **w,long *n_out,double sigma,double **err,ofstream &logout)
{
  double dof,res,res_vor;
  long *ind,ldum,k,l,m,o;
  double ***S,**SQ,**R;
  double *y,**gam,**x_out,**y_out,**g_out;

    
 
  ind=lvector(1,order);
  S=f3tensor(1,order,1,n,1,n);
  y=dvector(1,n);
  gam=dmatrix(1,order,1,n);
  x_out=dmatrix(1,order,1,n);
  y_out=dmatrix(1,order,1,n);
  g_out=dmatrix(1,order,1,n);

  R=dmatrix(1,n,1,n);
  SQ=dmatrix(1,n,1,n);

  for(o=1;o<=order;o++)
    {
      ind[o]=0;
    }
      
  for(k=1;k<=n;k++)
    {
      cerr << "Calculating error band(s) - point no. "  << k << "\r";
      //Einheitsvektor in y
      for(l=1;l<=n;l++)
	{
	  if(l==k)
	    y[l]=1.;
	  else
	    y[l]=0.;
	}

      res=10;
      //Backfitting
      for(l=1;l<=10000;l++)
	{	  
	  res_vor=res;
	  res=back_fit(x,y,n,order,alpha,g,gam,x_out,y_out,w,g_out,n_out,1);
	  
	  if(fabs(res_vor-res)<1e-3)
	    break;
	}

      //R-Matrix fuer Dof, wichtig : hier g statt g_out benutzen -> tied ranks
      for(l=1;l<=n;l++)
	{
	  R[l][k]=0.;
	  for(m=1;m<=order;m++)
	    R[l][k]+=g[m][l];
	}

      for(o=1;o<=order;o++)
	{
	  ldum=1;
	  for(l=k+1;l<=n;l++)
	    {
	      if(x[o][l]==x[o][k])
		ldum++;
	    }
	  if(ldum==1)
	    {
	      ind[o]++;
	      for(m=1;m<=n_out[o];m++)
		{
		  S[o][ind[o]][m]=g_out[o][m];	  
		}
	    }
	}
    }
  
  //Berechnung der Diagonalelemente der Kovarianzmatrix
  
  
  for(o=1;o<=order;o++)
    {
      for(k=1;k<=ind[o];k++)
	{
	  SQ[o][k]=0.;
	  for(l=1;l<=ind[o];l++)
	    {
	      SQ[o][k]+=S[o][l][k]*S[o][l][k];
	    }
	}
    }
 

  //Berechnung von dof
  dof=0.;
  for(k=1;k<=n;k++)
    {
      dof+=R[k][k];
    }
  sigma*=1./((double)n-dof);

  cerr << "Effective degrees of freedom :  " << dof << endl;
  cerr << "Estimated variance :  " << sqrt(sigma) << "\n\n\n";

  logout << endl;
  logout << "Effective degrees of freedom :  " << dof << endl;
  logout << "Estimated variance :  " << sqrt(sigma) << "\n\n";

  for(o=1;o<=order;o++)
    {
      for(k=1;k<=n_out[o];k++)
	{
	  err[o][k]=sqrt(sigma*SQ[o][k]);
	}
    }
  

  free_f3tensor(S,1,order,1,n,1,n);
  free_dvector(y,1,n);
  free_dmatrix(gam,1,order,1,n);
  free_dmatrix(x_out,1,order,1,n);
  free_dmatrix(y_out,1,order,1,n);
  free_dmatrix(g_out,1,order,1,n);
  
  
}
