#ifndef NR_H_TMP
#define NR_H_TMP
#include<fstream>

using namespace std;

// #####################################
// #                                   #
// #       Klasse f. Bandmatizen       #
// #                                   #
// #####################################


class Bandmatrix
{
private:
  long n;  //Anz. der Datenpunkte
  long bw; //Anz. der Nebendiagonalen
  double **band;
public:
  Bandmatrix(long ndata,long bandbreite); //Konstruktor
  double r(long i,long j); //Lesen
  void w(double value,long i,long j);  //Schreiben
  ~Bandmatrix();  //Destruktor
};

// #####################################
// #                                   #
// #      Splines Prototypes           #
// #                                   #
// #####################################


//Splines Berechnen
void splines(double *x,double *y,double *w,long n,double alpha,double *g,double *gam);
//Backfitting mittels Gauss-Seidel und Smoothing Splines als Glaetter
double back_fit(double **x,double *y,long n,long order,double *alpha,double **g,double **gam,double **x_out,double **y_out,double **w,double **g_out,long *n_out,int flag);
//Splines Plotten->interpolieren
double plot_splines(double *x,double *g,double *gam,long n,double t);
//Ableitung der Splines plotten
double deriv_splines(double *x,double *g,double *gam,long n,double t);
void error_band(double *x,double *y,long n,double *w,double alpha,double sigma,double *err);
double back_fit_kernel(double **x,double *y,long n,long order,double *h,double **g,double **x_out,double **y_out,double **w,double **g_out,long *n_out);
void kernel_smooth(double *x,double *y,long n,double h,double *w,double *g);
double K(double u,double h);
void splines_gcv(double *x,double *y,double *w,long n,double *alpha,double *g,double *gam);
long make_w(double *xs,double*ys,long *ar,long n,double *w);
void error_band(double **x,long n,long order,double *alpha,double **g,double **w,long *n_out,double sigma,double **err,ofstream &logout);
#endif

#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))

static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
        (dminarg1) : (dminarg2))

static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

static float minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))

static long lmaxarg1,lmaxarg2;
#define LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ?\
        (lmaxarg1) : (lmaxarg2))

static long lminarg1,lminarg2;
#define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?\
        (lminarg1) : (lminarg2))

static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void nrerror(string error_text);
int *ivector(long nl, long nh);
unsigned char *cvector(long nl, long nh);
long *lvector(long nl, long nh);
double *dvector(long nl, long nh);
char **cmatrix(long nrl, long nrh, long ncl, long nch);
long **lmatrix(long nrl, long nrh, long ncl, long nch);
float **matrix(long nrl, long nrh, long ncl, long nch);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);
float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl);
float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch);
double ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_ivector(int *v, long nl, long nh);
void free_cvector(unsigned char *v, long nl, long nh);
void free_lvector(long *v, long nl, long nh);
void free_dvector(double *v, long nl, long nh);
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_lmatrix(long **m, long nrl, long nrh, long ncl, long nch);
void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch);
void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch);
void free_f3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh);
#endif /* _NR_UTILS_H_ */

// #####################################
// #                                   #
// #        Prototypes for RECIPES     #
// #                                   #
// #####################################

#ifndef _NR_H_
#define _NR_H_

#ifndef _FCOMPLEX_DECLARE_T_
typedef struct FCOMPLEX {float r,i;} fcomplex;
#define _FCOMPLEX_DECLARE_T_
#endif /* _FCOMPLEX_DECLARE_T_ */

#ifndef _ARITHCODE_DECLARE_T_
typedef struct {
	unsigned long *ilob,*iupb,*ncumfq,jdif,nc,minint,nch,ncum,nrad;
} arithcode;
#define _ARITHCODE_DECLARE_T_
#endif /* _ARITHCODE_DECLARE_T_ */

#ifndef _HUFFCODE_DECLARE_T_
typedef struct {
	unsigned long *icod,*ncod,*left,*right,nch,nodemax;
} huffcode;
#define _HUFFCODE_DECLARE_T_
#endif /* _HUFFCODE_DECLARE_T_ */

#include <stdio.h>

double gammln(double xx);
double gammp(double a, double x);
double gammq(double a, double x);
void gcf(double *gammcf, double a, double x, double *gln);
void gser(double *gamser, double a, double x, double *gln);
void indexx(unsigned long n, double arr[], long indx[]);
double ran1(long *idum);
void choldc(double **a, int n);
void Choldc(Bandmatrix *a, int n, double p[]);
void cholsl(double **a, int n, double p[], double b[], double x[]);
void Cholsl(Bandmatrix *a, int n, double p[], double b[], double x[]);

#endif /* _NR_H_ */
