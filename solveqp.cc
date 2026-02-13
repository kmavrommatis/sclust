//########################################################################################
//
//Sclust: Mutation Clustering using Smoothing Splines
//Copyright (C) 2017, Peifer Lab, University of Cologne
//
//This program is free software: you can redistribute it and/or modify it under 
//the terms of the GNU General Public License as published by the Free Software 
//Foundation, either version 3 of the License, or (at your option) any later version.
//
//This program is distributed in the hope that it will be useful, but 
//WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
//FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License along with this
//program. If not, see <http://www.gnu.org/licenses/>.
//########################################################################################


#include "solveqp.h"
#include <iomanip>

using namespace std;

//macro for 2d-access on 1d-array, rs is rowsize, x is column, y is row
  #define at(x,y,rs) (((x)*(rs)+(y)))

void clear_dvector(double *v,long n)
{
  for(long i=1;i<=n;i++)
    v[i]=0;
}

void clear_dmatrix(double **m,long nr,long nc)
{
  for(long i=1;i<=nr;i++)
    {
      for(long j=1;j<=nc;j++)
	{
	  m[i][j]=0;
	}
    }
}

double rp(double val,long prec)
{
  if(fabs(val) <= pow(10,-prec))
    return(0);
  else
    return(val);
}
#ifdef LAPACK
void choldc_inv(long n,double **G,double **Ginv)
{

   //copy 2d-array into 1d-array
  prec_t *G_p, *Ginv_p;
  G_p = (double *)mkl_malloc( n*n*sizeof( double ), 64 );

  for (size_t i=0;i<n;++i)
    {
      for (size_t j=0;j<n;++j)
	{
	  G_p[at(i,j,n)]=G[i+1][j+1];
	}
    }

  //   std::cout << "Using Lapack for Inversion\n";
  //  might use dgetri in the future because of band properties
  lapack_int tmp=LAPACKE_dpotrf(LAPACK_ROW_MAJOR,'U',(lapack_int) n,G_p,(lapack_int) n);      //choles
  if (tmp!=0) std::cout << "dpotrf: faulty!\n";
  lapack_int tmp2=LAPACKE_dpotri(LAPACK_ROW_MAJOR,'U',(lapack_int) n,G_p,(lapack_int) n);  //Compute i
  if (tmp2!=0) std::cout << "dpotri: faulty!\n";
  //Only upper triangular matrices is inverse - copy to lower trian.
#pragma GCC ivdep
#pragma ivdep
  for (size_t r=0;r<n;++r)
    {
      for (size_t c=0;c<n;++c) 
	{
	  if (r<=c)
            Ginv[r+1][c+1]=G_p[at(r,c,n)];
	  else
            Ginv[r+1][c+1]=G_p[at(c,r,n)];
	}
    }

  mkl_free(G_p);
}

#else
void choldc_inv(long n,double **G,double **Ginv)
{
  double *p;
  double **Gtmp;
  double sum;
  long i,j,k;

  p=dvector(1,n);
  Gtmp=dmatrix(1,n,1,n);
  

  for(i=1;i<=n;i++)
    {
      for(j=1;j<=n;j++)
	Gtmp[i][j]=G[i][j];
    }

  choldc(Gtmp,n);

  for(i=1;i<=n;i++)
    p[i]=Gtmp[i][i];
  
  for (i=1;i<=n;i++) 
    { 
      Gtmp[i][i]=1.0/p[i];
      for (j=i+1;j<=n;j++) 
	{
	  sum=0.0;
	  for (k=i;k<j;k++) 
	    sum -= Gtmp[j][k]*Gtmp[k][i];
	  Gtmp[j][i]=sum/p[j];
	} 
    }
  
   clear_dmatrix(Ginv,n,n);
   for(i=1;i<=n;i++)
    {
      for(j=1;j<=n;j++)
	{
	  for(k=1;k<=n;k++)
	    Ginv[i][j]+=Gtmp[k][i]*Gtmp[k][j];
	}
    }
  
  free_dvector(p,1,n);
  free_dmatrix(Gtmp,1,n,1,n);
}
#endif   //LAPACK

void initmat(prec_t *MAT, size_t r, size_t c, prec_t val) {
#pragma GCC ivdep
#pragma ivdep
   for (size_t i = 0; i < r * c; ++i)
      MAT[i] = val;
}

void outmat(prec_t *M, size_t rows, size_t cols, char *descr) {
  std::cout << "1D: " << descr << ":" << std::endl;
  std::cout << std::setprecision(2);
   for (size_t r=0;r<rows; ++r) {
     for (size_t c=0;c<cols; ++c) {
         std::cout << M[at(r,c,cols)] << ' ';
     }
     std::cout << std::endl;
   }
}

void outmat_2d(prec_t **M, size_t rows, size_t cols, char *descr) {
  std::cout << "2D: " << descr << ":" << std::endl;
  std::cout << std::setprecision(2);
  for (size_t r=1;r<=rows; ++r) {
      for (size_t c=1;c<=cols; ++c) {
         std::cout << M[r][c] << ' ';
      }
      std::cout << std::endl;
   }
}

void outvec_2d(prec_t *M, size_t rows, char *descr) {
    std::cout << "2D: " << descr << ":" << std::endl;
    std::cout << std::setprecision(2);
    for (size_t r=1;r<=rows; ++r) {
      std::cout << M[r] << std::endl;
    }
//    std::cout << std::endl;
      }

#ifdef LAPACK

void solve_sub(long n, double *Ginv_M, double *a, vector<long> &act_v, size_t act_v_size, double *u_v, double *S_M)
{
   size_t na;

   initmat(S_M,n,n,0);  //really necessary?
   initmat(u_v,n,1,0);//really necessary?

   if(act_v_size==0)//solve the unconstrained problem
   {
#pragma GCC ivdep
#pragma ivdep
      for (size_t i=0;i<n*n;++i)
         S_M[i]=-Ginv_M[i];
   }
   else
   {
      na=act_v_size;

      //create projection matrix N
      prec_t *M,*TMP_na_n,*Minv,*N,*Nstar;

      N= (double *)mkl_malloc( n*na*sizeof( double ), 64 );
      initmat(N,n,na,0);

      Minv= (double *)mkl_malloc( na*na*sizeof( double ), 64 );

      TMP_na_n= (double *)mkl_malloc( na*n*sizeof( double ), 64 );

      Nstar= (double *)mkl_malloc( na*n*sizeof( double ), 64 );
      initmat(Nstar,na,n,0);

#pragma GCC ivdep
#pragma ivdep
      for(size_t j=0;j<na;++j)
         N[at(act_v[j],j,na)]=1;    //careful, still correct?

      // M=N^T%*%Ginv%*%N     //AT-OK
      //First: TMP_na_n=N%*%Ginv
      matmul('T','N',N,n,na,Ginv_M,n,n,TMP_na_n,na,n);
      //Last: M=TMP_na_n%*%N
      matmul('N','N',TMP_na_n,na,n,N,n,na,Minv,na,na);

      lapack_int tmp=LAPACKE_dpotrf(LAPACK_ROW_MAJOR,'U',(lapack_int) (na),Minv,(lapack_int)(na));//choleski-Fak.: for symmetric, hermitian
      if (tmp!=0) std::cout << "dpotrf: faulty!\n";
      lapack_int tmp2=LAPACKE_dpotri(LAPACK_ROW_MAJOR,'U',(lapack_int) (na),Minv,(lapack_int)(na));//Compute inverse of sem, herm., positive-def. matrix using Cholesky fac.
      if (tmp2!=0) std::cout << "dpotri: faulty!\n";
      //Only upper triangular matrices is inverse - copy to lower trian.
#pragma GCC ivdep
#pragma ivdep
      for (size_t r=0;r<na;++r) {
         for (size_t c=r+1;c<na;++c)
         Minv[at(c,r,na)]=Minv[at(r,c,na)];
      }

      //Nstar = Minv %*% N^T %*% Ginv     //AT-OK
      //first: TMP=Minv %*% N
      matmul('N','T',Minv,na,na,N,n,na,TMP_na_n,na,n);
      //Last: Nstar = TMP_na_na %*% Ginv
      matmul('N','N',TMP_na_n,na,n,Ginv_M,n,n,Nstar,na,n);

      //compute Lagrange multipliers      //maybe BLAS usable?
#pragma GCC ivdep
#pragma ivdep
      for(size_t i=0;i<na;++i) {
         for(size_t j=0;j<n;++j) {
            u_v[act_v[i]]+=Nstar[at(i,j,n)]*a[j];
         }
      }
      //compute solution matrix S
      M= (double *)mkl_malloc( n*n*sizeof( double ), 64 );
      initmat(M,n,n,0);

      //M = -1 * N %*% Nstar     //AT-OK
      matmul('N','N',N,n,na,Nstar,na,n,-1,0,M,n,n);

#pragma GCC ivdep
#pragma ivdep
      for(size_t i=0;i<n;++i)
         M[at(i,i,n)]+=1;

      //S = -1 * Ginv %*% M      //AT-OK
      matmul('N','N',Ginv_M,n,n,M,n,n,-1,0,S_M,n,n);

      mkl_free(TMP_na_n);
      mkl_free(N);
      mkl_free(M);
      mkl_free(Minv);
      mkl_free(Nstar);
   }
}

#else    //Means no LAPACK..
void solve_sub(long n,double **Ginv,double *a,vector<long> act,double *u,double **S)
{
  long i,j,k,l,m;
  long na;
  double **M,**Minv,**N,**Nstar;

  clear_dmatrix(S,n,n);
  clear_dvector(u,n);

  if(act.size()==0) //solve the unconstrained problem
    {
      for(i=1;i<=n;i++)
	{
	  for(j=1;j<=n;j++)
	    {
	      S[i][j]=-Ginv[i][j];
	    }
	}
    }
  else
    {
      na=act.size();

      //create projection matrix N
      N=dmatrix(1,n,1,na);
      M=dmatrix(1,na,1,na);
      Minv=dmatrix(1,na,1,na);
      Nstar=dmatrix(1,na,1,n);

      clear_dmatrix(N,n,na);
      clear_dmatrix(M,na,na);
      clear_dmatrix(Nstar,na,n);

      for(j=1;j<=na;j++)
	N[act[j-1]][j]=1;
       
      for(i=1;i<=na;i++)
	{
	 for(j=1;j<=na;j++)
	    {
	      for(k=1;k<=n;k++)
		{
		  for(l=1;l<=n;l++)
		    {
		      M[i][j]+=N[l][i]*Ginv[l][k]*N[k][j];
		    }
		}
	    }
	}

      choldc_inv(na,M,Minv);
      free_dmatrix(M,1,na,1,na);

      for(i=1;i<=na;i++)
	{
	  for(j=1;j<=n;j++)
	    {
	      for(k=1;k<=n;k++)
		{
		  for(l=1;l<=na;l++)
		    {
		      Nstar[i][j]+=Minv[i][l]*N[k][l]*Ginv[k][j];
		    }
		}
	    }
	}

      //compute Lagrange multipliers
      for(i=1;i<=na;i++)
	{
	  for(j=1;j<=n;j++)
	    {
	      u[act[i-1]]+=Nstar[i][j]*a[j];
	    }
	}

      //compute solution matrix S
      M=dmatrix(1,n,1,n);
      clear_dmatrix(M,n,n);
      for(i=1;i<=n;i++)
	{
	  for(j=1;j<=n;j++)
	    {
	      for(k=1;k<=na;k++)
		{
		  M[i][j]+=-N[i][k]*Nstar[k][j];
		}
	    }
	}

      for(i=1;i<=n;i++)
	M[i][i]+=1;

      for(i=1;i<=n;i++)
	{
	  for(j=1;j<=n;j++)
	    {
	      for(k=1;k<=n;k++)
		{
		  S[i][j]+=-Ginv[i][k]*M[k][j];
		}
	    }
	}

      free_dmatrix(N,1,n,1,na);
      free_dmatrix(Nstar,1,na,1,n);
      free_dmatrix(M,1,n,1,n);
      free_dmatrix(Minv,1,na,1,na);
    }
}
#endif      //LAPACK


#ifdef LAPACK

void cpymat(double *SOURCE, double *TARGET, size_t row, size_t col) {
     cblas_dcopy(row*col,SOURCE,1,TARGET,1);
}

void solveqp(long n,double **G,double *a,double **S,long maxit) {
//wrapper for 2d-Arrays
   prec_t *G_Mat, *a_vec, *S_Mat;
   G_Mat= (double *)mkl_malloc( n*n*sizeof( double ), 64 );
   a_vec= (double *)mkl_malloc( n*sizeof( double ), 64 );
   S_Mat= (double *)mkl_malloc( n*n*sizeof( double ), 64 );

   for (size_t i=0;i<n;++i)
      for (size_t j=0;j<n;++j)
         G_Mat[at(i,j,n)]=G[i+1][j+1];

   for (size_t i=0;i<n;++i)
      a_vec[i]=a[i+1];

   solveqp_1D(n, G_Mat, a_vec, S_Mat, maxit);

   for (size_t i=0;i<n;++i)
      for (size_t j=0;j<n;++j)
         S[i+1][j+1]=S_Mat[at(i,j,n)];

   mkl_free(G_Mat);
   mkl_free(a_vec);
   mkl_free(S_Mat);
}

void solveqp_1D(long n,double *G,double *a,double *S,long maxit)
{  
   //solves the problem min(a^T x+0.5 x^T G x) subject to x>=0
   long prec=6;
   long nx;
   vector<long> act;//active constraints
   
   long it=0;

   prec_t *Ginv,*u,*x,*x_pre;

   Ginv= (double *)mkl_malloc( n*n*sizeof( double ), 64 );
   initmat(Ginv,n,n,0);
   u= (double *)mkl_malloc( n*sizeof( double ), 64 ); //Lagrange Multipliers
   initmat(u,n,1,0);    //is first used in ..._sub and initialized there so probalby no need to init here
   x= (double *)mkl_malloc( n*sizeof( double ), 64 ); //Solution Vector
   initmat(x,n,1,0);
   x_pre= (double *)mkl_malloc( n*sizeof( double ), 64 );
   initmat(x_pre,n,1,0);
   
   cpymat(G,Ginv,n,n);
   
   lapack_int tmp=LAPACKE_dpotrf(LAPACK_ROW_MAJOR,'U',(lapack_int) (n),Ginv,(lapack_int)(n));//choleski-Fak.: for symmetric, hermitian
   if (tmp!=0) std::cout << "dpotrf: faulty!\n";
   lapack_int tmp2=LAPACKE_dpotri(LAPACK_ROW_MAJOR,'U',(lapack_int) (n),Ginv,(lapack_int)(n));//Compute inverse of sem, herm., positive-def. matrix using Cholesky fac.
   if (tmp2!=0) std::cout << "dpotri: faulty!\n";
   //Only upper triangular matrices is inverse - copy to lower trian.
#pragma GCC ivdep
#pragma ivdep
   for (size_t r=0;r<n;++r) {
      for (size_t c=r+1;c<n;++c)
         Ginv[at(c,r,n)]=Ginv[at(r,c,n)];
   }

   //initialization
   act.clear();
   solve_sub(n,Ginv,a,act,act.size(),u,S);//careful: act might be not be acceptable as vector
   initmat(x,n,1,0);
   initmat(x_pre,n,1,0);
   
   //x = S %*% a
   matmul('N','N',S,n,n,a,n,1,x,n,1);
   
   //iterations
   nx=n;
   while(nx!=0) {
      it++;
   
      cpymat(x,x_pre,n,1);

      //find set of active constraints to test
      act.clear(); nx=0;
      for(size_t i=0;i<n;++i) {
         if(rp(x[i],prec) < 0.0)
            nx++;
         if(rp(x[i],prec) < 0.0 || rp(u[i],prec) > 0.0)
            act.push_back(i);
      }  
            
      //...and solve again
      solve_sub(n,Ginv,a,act,act.size(),u,S);
      initmat(x,n,1,0);
      
      //x = S %*% a
      matmul('N','N',S,n,n,a,n,1,x,n,1);
      
      if(it==maxit) {
         cerr << "Error: QP has reached maximal number of iterations currently set to: " << maxit << ".\n";
         exit(2);
      }
   }

   mkl_free(Ginv);
   mkl_free(u);
   mkl_free(x);
   mkl_free(x_pre);
}
#else
void solveqp(long n,double **G,double *a,double **S,long maxit)
{
  //solves the problem min(a^T x+0.5 x^T G x) subject to x>=0
  double **Ginv;
  long prec=6;
  double *u; //Lagrange multipliers
  double *x,*x_pre; //solution

  long nx;
  vector<long> act; //active constraints

  long i,j,it=0;

  Ginv=dmatrix(1,n,1,n);  
  u=dvector(1,n);
  x=dvector(1,n);
  x_pre=dvector(1,n);

  choldc_inv(n,G,Ginv);
  
  //initialization
  act.clear();
  solve_sub(n,Ginv,a,act,u,S);
  clear_dvector(x_pre,n);
  clear_dvector(x,n);
  for(i=1;i<=n;i++)
    {
      for(j=1;j<=n;j++)
	{
	  x[i]+=S[i][j]*a[j];
	}
    }

  //iterations
  nx=n;
  while(nx!=0)
    {
      it++;
      for(i=1;i<=n;i++)
	x_pre[i]=x[i];

      //find set of active constraints to test
      act.clear(); nx=0;
      for(i=1;i<=n;i++)
	{
	  if(rp(x[i],prec) < 0.0)
	    nx++;

	  if(rp(x[i],prec) < 0.0 || rp(u[i],prec) > 0.0)
	    {
	      act.push_back(i);
	    }
	}

      //...and solve again
      solve_sub(n,Ginv,a,act,u,S);
      clear_dvector(x,n);
      for(i=1;i<=n;i++)
	{
	  for(j=1;j<=n;j++)
	    {
	      x[i]+=S[i][j]*a[j];
	    }
	}

      if(it==maxit)
	{
	  cerr << "Error: QP has reached maximal number of iterations currently set to: " << maxit << ".\n";
	  exit(2);
	}
    }

  free_dmatrix(Ginv,1,n,1,n);
  free_dvector(u,1,n);
  free_dvector(x,1,n);
  free_dvector(x_pre,1,n);
}
#endif      //LAPACK
