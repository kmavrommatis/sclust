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

#ifndef QP_H
#define QP_H

#include <iostream>
#include <limits>    //still necessary?
#include <iomanip>   //included for setpreision
#include "math.h"
#include <stdio.h>
#include "solveqp.h"

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

typedef double prec_t;

void vinit(prec_t *M,size_t n);
inline prec_t pnorm(prec_t q, prec_t mean, prec_t sd);
inline prec_t calc_a (size_t i, size_t j, prec_t *x, prec_t sig);
void diag(prec_t *MAT, prec_t val, size_t n);
void spline_deconv_qp (prec_t *Dmat_p, prec_t *dvec_p, prec_t *S_p, prec_t *x_p, size_t n);
void matmul(char tA, char tB, prec_t *A, size_t A_rows, size_t A_cols, prec_t *B, size_t B_rows, size_t B_cols, prec_t *C, size_t C_rows, size_t C_cols);
void matmul(char tA, char tB, prec_t *A, size_t A_rows, size_t A_cols, prec_t *B, size_t B_rows, size_t B_cols, prec_t Sca_Mul, prec_t Sca_Add, prec_t *C, size_t C_rows, size_t C_cols);
void spline_deconv(prec_t *x, size_t x_s, prec_t *y, size_t y_s, prec_t *w, prec_t *g, prec_t *gam_p, prec_t alpha);
#endif
