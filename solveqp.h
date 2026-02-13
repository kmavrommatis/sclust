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

//set some include guards:
#ifndef SOLVEQP_H
#define SOLVEQP_H

#include <stdio.h>
#include <stdint.h>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <getopt.h>
#include <math.h>
#include <time.h>
#include <map>
#include <list>
#include <sstream>
#include "qp++.h"

#include "numerics/nr.h"

//#define MKL     uncomment or define MKL at compile time, to use MKL instead of openBLAS
#ifdef MKL
#include "mkl.h"
#else
#include "openBLAS/include/cblas.h"
#ifdef LAPACK
#include "openBLAS/include/lapacke.h"
#endif
#define mkl_malloc(x,y) malloc( (x) )
#define mkl_free(x) free( (x) )
#define lapack_int int
#endif

//macro for 2d-access on 1d-array, rs is rowsize, x is column, y is row
#define at(x,y,rs) (((x)*(rs)+(y)))

typedef double prec_t;


void clear_dvector(double *v,long n);
void clear_dmatrix(double **m,long nr,long nc);
double rp(double val,long prec);
void choldc_inv(long n,double **G,double **Ginv);
double l1_norm(long n,double *x1,double *x2);
void initmat(prec_t *MAT, size_t r, size_t c, prec_t val);
void solveqp(long n,double **G,double *a,double **S,long maxit);
#ifdef LAPACK
void cpymat(double *SOURCE, double *TARGET, size_t row, size_t col);
void solveqp_1D(long n,double *G,double *a,double *S,long maxit);
void solve_sub(long n, double *Ginv_M, double *a, vector<long> &act_v, size_t act_v_size, double *u_v, double *S_M);
#else
void solve_sub(long n,double **Ginv,double *a,vector<long> act,double *u,double **S);
#endif   //LAPACK
#endif   //Header-Guard

