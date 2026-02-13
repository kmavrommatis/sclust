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

#include "qp++.h"

//#define DEBUG		//Uncomment for production code
#ifdef DEBUG
#define OUTPUT		//if DEBUG, uncomment this to prevent Outputs of matrices
#endif
#define OUTPUT_PREC 2

#define min(x,y) (((x) < (y)) ? (x) : (y))
//macro for 2d-access on 1d-array, rs is rowsize, x is column, y is row
#define at(x,y,rs) (((x)*(rs)+(y)))

inline prec_t pnorm(prec_t q, prec_t mean, prec_t sd) 
{
  prec_t tmp=q/(sqrt((double) 2)*sd);
  return ((prec_t) 0.5 * (prec_t) erf((double) tmp))+(prec_t)0.5;
}

inline prec_t calc_a (size_t i, size_t j, prec_t *x, prec_t sig) 
{
  return (pnorm(x[j]-x[i],0,sig)-pnorm(x[j]-x[i+1],0,sig));
}

void diag(prec_t *MAT, prec_t val, size_t n) 
{
  // init all values with zero
  for (size_t i=0;i<n*n;++i) {
    MAT[i]=0;
  }
  // set diag to val
  for (size_t i=0;i<n;++i) 
    {      
      MAT[at(i,i,n)]=val;
    }
}

void vinit(prec_t *M,size_t n)
{
#pragma GCC ivdep
#pragma ivdep
  for (size_t i=0;i<n;i++)
    M[i]=0;
}

#ifdef LAPACK
void spline_deconv_qp (prec_t *Dmat_p, prec_t *dvec_p, prec_t *S_p, prec_t *x_p, size_t n, long max_qp_iter) {
   long prec=6;

   prec_t *a= (double *)mkl_malloc( n*sizeof( double ), 64 );;
//   initmat(S_p,n,n,0);   //AT: Not necessary, since initialized in solve_sub(...)
//   initmat(x_p,n,1,0);  //AT: Not necessary, since overwirtten in matmul(...)

#pragma GCC ivdep
#pragma ivdep
   for(size_t i=0;i<n;++i)
      a[i]=(-1)*dvec_p[i];


   solveqp_1D(n,Dmat_p,a,S_p,max_qp_iter);

   matmul('N','N',S_p,n,n,a,n,1,x_p,n,1);

   mkl_free(a);
}

#else
void spline_deconv_qp (prec_t *Dmat_p, prec_t *dvec_p, prec_t *S_p, prec_t *x_p, size_t n, long max_qp_iter) {

  long prec=6;

  prec_t **G,*a,*x,**S;

  G=dmatrix(1,n,1,n);
  S=dmatrix(1,n,1,n);
  a=dvector(1,n);
  x=dvector(1,n);
	
#pragma GCC ivdep
  for(size_t i=1;i<=n;++i)
    a[i]=(-1)*dvec_p[i-1];

#pragma GCC ivdep
  for(size_t i=1;i<=n;++i)
    {
      for(size_t j=1;j<=n;++j)
	{
	  G[i][j]=Dmat_p[at((i-1),(j-1),n)];
	}
    }
  
  clear_dmatrix(S,n,n);
  
  //solves the problem min(a^T x+0.5 x^T G x) subject to x>=0
  solveqp(n,G,a,S,max_qp_iter);
  
  clear_dvector(x,n);
  for(size_t i=1;i<=n;++i)
    {
      for(size_t j=1;j<=n;++j)
	{
	  x[i]+=S[i][j]*a[j];
	}
    }
  
  //cpy back into one-dimensional layout: S
#pragma GCC ivdep
  for (size_t i=1;i<=n;++i) {
    for (size_t j=1;j<=n;++j) {
      S_p[at((i-1),(j-1),n)]=S[i][j];
    }
  }
  //cpy back into one-dimensional layout: x
#pragma GCC ivdep
  for (size_t i=1;i<=n;++i) {
    x_p[i-1]=x[i];
  }
  
  free_dmatrix(G,1,n,1,n);
  free_dmatrix(S,1,n,1,n);
  free_dvector(a,1,n);
  free_dvector(x,1,n);
}
#endif      //LAPACK

void matmul(char tA, char tB, prec_t *A, size_t A_rows, size_t A_cols, prec_t *B, size_t B_rows, size_t B_cols, prec_t *C, size_t C_rows, size_t C_cols) {
   matmul(tA, tB, A, A_rows, A_cols, B, B_rows, B_cols, 1, 0, C, C_rows, C_cols);
}

void matmul(char tA, char tB, prec_t *A, size_t A_rows, size_t A_cols, prec_t *B, size_t B_rows, size_t B_cols, prec_t Sca_Mul, prec_t Sca_Add, prec_t *C, size_t C_rows, size_t C_cols) {
  
  // Form appropriate dgemm Statements for matmul function call; Always first check, if call is defined 
  if (tA=='N' && tB=='N') {
    if (A_cols==B_rows && A_rows==C_rows && B_cols==C_cols) {
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, A_rows, B_cols, A_cols, Sca_Mul, A, A_cols, B, B_cols, Sca_Add, C, C_cols);
    }
    else {
      std::cerr << "MatMul(N,N): For these Matrices, matrix multiplication is not defined! This is a coding Error. Exiting..\n";
      exit(1);
    }
  }
  else if (tA=='N' && tB=='T') {
    if (A_cols==B_cols && A_rows==C_rows && B_rows==C_cols) {
      cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, A_rows, B_rows, A_cols, Sca_Mul, A, A_cols, B, B_cols, Sca_Add, C, C_cols);
    }
    else {
      std::cerr << "MatMul(N,T): For these Matrices, matrix multiplication is not defined! This is a coding Error. Exiting..\n";
      exit(1);
    }
  }
  else if (tA=='T' && tB=='N') {
    if (A_rows==B_rows && A_cols==C_rows && B_cols==C_cols) {
      cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, A_cols, B_cols, A_rows, Sca_Mul, A, A_cols, B, B_cols, Sca_Add, C, C_cols);
    }
    else {
      std::cerr << "MatMul(T,N): For these Matrices, matrix multiplication is not defined! This is a coding Error. Exiting..\n";
      exit(1);
    }
  }
  else if (tA=='T' && tB=='T') {
    if (A_rows==B_cols && A_cols==C_rows && B_rows==C_cols) {
      cblas_dgemm(CblasRowMajor, CblasTrans, CblasTrans, A_cols, B_rows, A_rows, Sca_Mul, A, A_cols, B, B_cols, Sca_Add, C, C_cols);
    }
    else {
      std::cerr << "MatMul(T,T): For these Matrices, matrix multiplication is not defined! This is a coding Error. Exiting..\n";
      exit(1);
    }
  }
  else {
    std::cerr << "MatMul: Malfromed Function call! This is a coding Error. Exiting..\n";
    exit(1);
  }
}

void spline_deconv(prec_t *x, size_t x_s, prec_t *y, size_t y_s, prec_t *w, prec_t *g,prec_t *GAM, prec_t alpha, long max_qp_iter) {
  size_t n=x_s;
  double gcv,trS;
  prec_t *Q,*R,*W,*TMP,*C,*A,*A_TMP,*Dmat,*Amat,*dvec,*RTST,*RBAK,*GAM_TMP,*GAM_SUB_TMP,*S_tmp;
  //Allocating memory for matrices aligned on 64-byte boundary for better Performance
  Q = (double *)mkl_malloc( n*(n-2)*sizeof( double ), 64 );
  R = (double *)mkl_malloc( (n-2)*(n-2)*sizeof( double ), 64 );
  TMP = (double *)mkl_malloc( (n)*(n-2)*sizeof( double ), 64 );
  C = (double *)mkl_malloc( (n)*(n)*sizeof( double ), 64 );
  A = (double *)mkl_malloc( n*n*sizeof( double ), 64 );
  A_TMP = (double *)mkl_malloc( n*n*sizeof( double ), 64 );
  Dmat = (double *)mkl_malloc( n*n*sizeof( double ), 64 );
  dvec = (double *)mkl_malloc( n*sizeof( double ), 64 );
  S_tmp = (double *)mkl_malloc( n*n*sizeof( double ), 64 );
  GAM_SUB_TMP = (double *)mkl_malloc( (n-2)*(n)*sizeof( double ), 64 );
  GAM_TMP = (double *)mkl_malloc( (n-2)*sizeof( double ), 64 );
  
  if (Q == NULL || R == NULL || C == NULL || A == NULL || TMP == NULL || A_TMP == NULL || Dmat == NULL || dvec == NULL || GAM_SUB_TMP == NULL || GAM_TMP == NULL) {
    std::cerr << "\n ERROR: Can't allocate memory for matrices. Aborting... \n\n";
    mkl_free(Q);
    mkl_free(R);
    mkl_free(C);
    mkl_free(A);
    mkl_free(TMP);
    mkl_free(A_TMP);
    mkl_free(Dmat);
    mkl_free(dvec);
    mkl_free(GAM_SUB_TMP);
    mkl_free(GAM_TMP);
    exit(1);
  }
  
  vinit(Q,n*(n-2));
  vinit(R,(n-2)*(n-2));

#pragma GCC ivdep
  for (size_t i=0;i<n-2;i++) 
    {
      Q[at(i,i,n-2)]=1.0/(x[i+1]-x[i]);
      Q[at(i+1,i,n-2)]=-1.0/(x[i+1]-x[i])-1.0/(x[i+2]-x[i+1]);
      Q[at(i+2,i,n-2)]=1.0/(x[i+2]-x[i+1]);
    }

#pragma GCC ivdep
#pragma ivdep
  for (size_t i=0;i<n-2;++i) 
    {
      R[at(i,i,n-2)]=((prec_t)1.0/(prec_t)3.0)*(x[i+2]-x[i]);
      if (i<n-3) 
	{
	  R[at(i+1,i,n-2)]=((prec_t)1.0/(prec_t)6.0)*(x[i+2]-x[i+1]);
	  R[at(i,i+1,n-2)]=((prec_t)1.0/(prec_t)6.0)*(x[i+2]-x[i+1]);
	}
    }

  //solve(R)-> R^-1
#ifdef LAPACK
  //   std::cout << "Using Lapack for Inversion\n";
  //	might use dgetri in the future because of band properties
  lapack_int tmp=LAPACKE_dpotrf(LAPACK_ROW_MAJOR,'U',(lapack_int) (n-2),R,(lapack_int)(n-2));		//choleski-Fak.: for symmetric, hermitian
  if (tmp!=0) std::cout << "dpotrf: faulty!\n";
  lapack_int tmp2=LAPACKE_dpotri(LAPACK_ROW_MAJOR,'U',(lapack_int) (n-2),R,(lapack_int)(n-2));	//Compute inverse of sem, herm., positive-def. matrix using Cholesky fac.
  if (tmp2!=0) std::cout << "dpotri: faulty!\n";
  //Only upper triangular matrices is inverse - copy to lower trian.
#pragma GCC ivdep
#pragma ivdep
  for (size_t r=0;r<n-2;++r) 
    {
      for (size_t c=r+1;c<n-2;++c)
	R[at(c,r,n-2)]=R[at(r,c,n-2)];
    }
#else
  //   std::cout << "Using choldc_inv for Inversion\n";
  double **R_t,**R_inv;
  R_t=dmatrix(1,n-2,1,n-2);
  R_inv=dmatrix(1,n-2,1,n-2);
  for (size_t r=0;r<n-2;r++) 
    {
      for (size_t c=0;c<n-2;c++)
	R_t[r+1][c+1]=R[at(r,c,n-2)];
    }
  choldc_inv(n-2,R_t,R_inv);
  for (size_t r=0;r<n-2;r++) 
    {
      for (size_t c=0;c<n-2;c++)
	R[at(r,c,n-2)]=R_inv[r+1][c+1];
    }
  free_dmatrix(R_t,1,n-2,1,n-2);
  free_dmatrix(R_inv,1,n-2,1,n-2);
#endif
  
  
  //C=Q%*%R^-1%*%t(Q)
  matmul('N','N',Q,n,n-2,R,n-2,n-2,TMP,n,n-2);
  matmul('N','T',TMP,n,n-2,Q,n,n-2,C,n,n);
  //maybe cblas_dsymm possible (if symmetric); maybe cblas_dhemm if hermitian
  
  prec_t sig;
#pragma GCC ivdep
#pragma ivdep
  for (size_t i=1;i<n-1;++i) 
    {
      for (size_t j=0; j<n;++j) 
	{
	  sig=w[i];
	  A[at(i,j,n)]=(prec_t)0.5*(calc_a(i,j,x,sig)+calc_a(i-1,j,x,sig));
	}
    }
#pragma GCC ivdep
#pragma ivdep
  for (size_t j=0;j<n;++j) 
    {
      A[at(0,j,n)]=(prec_t)0.5*(calc_a(0,j,x,sig));
      A[at(n-1,j,n)]=(prec_t)0.5*(calc_a(n-2,j,x,sig));
    }	
  
  //Dmat=A%*%t(A)+alpha*K; A_TMP=a%*%t(A)
  matmul('N','T',A,n,n,A,n,n,A_TMP,n,n);

  //Dmat=A_TMP+alpha*C
#pragma GCC ivdep
#pragma ivdep
  for (size_t i=0;i<n*n;++i) 
    {
      Dmat[i]=A_TMP[i]+alpha*C[i];
    }

  //dvec=t(y)%*%t(A)
  matmul('T','T',y,n,1,A,n,n,dvec,1,n);	
  
  //call quadratic solver; Solution into g vector and S_tmp Matrix
  spline_deconv_qp(Dmat, dvec, S_tmp, g, n, max_qp_iter);
  
  // gam=c(0,solve(R)%*%t(Q)%*%g,0)	with R still beeing R^-1=solve(R)
  // first: GAM_SUB_TMP=solve(R)%*%t(Q)
  matmul('N','T',R,n-2,n-2,Q,n,n-2,GAM_SUB_TMP,n-2,n);
  // then: GAM_TMP=GAM_SUB_TMP%*%g
  matmul('N','N',GAM_SUB_TMP,n-2,n,g,n,1,GAM_TMP,n-2,1);

  //GAM=c(0,GAM_TMP,0)
  GAM[0]=0;
  GAM[n-1]=0;
  for(size_t i=1;i<n-1;++i)
    GAM[i]=GAM_TMP[i-1];

  //clean up memory
  mkl_free(Q);
  mkl_free(R);
  mkl_free(C);
  mkl_free(A);
  mkl_free(TMP);
  mkl_free(A_TMP);
  mkl_free(Dmat);
  mkl_free(dvec);
  mkl_free(S_tmp);
  mkl_free(GAM_TMP);
  mkl_free(GAM_SUB_TMP);
}
