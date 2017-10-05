/*
 Genome-wide Efficient Mixed Model Association (GEMMA)
 Copyright (C) 2011  Xiang Zhou and Matthew Stephens
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>
#include <bitset>
#include <vector>
#include <map>
#include <set>
#include <cstring>
#include <cmath>
#include <stdio.h>
#include <stdlib.h> 

#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_cdf.h"

#ifdef FORCE_FLOAT
#include "mathfunc_float.h"
#else
#include "mathfunc.h"
#endif


using namespace std;



//calculate variance of a vector
double VectorVar (const gsl_vector *v)
{
	double d, m=0.0, m2=0.0;
	for (size_t i=0; i<v->size; ++i) {
		d=gsl_vector_get (v, i);
		m+=d;
		m2+=d*d;
	}
	m/=(double)v->size;
	m2/=(double)v->size;
	return m2-m*m;
}



//center the matrix G	
void CenterMatrix (gsl_matrix *G)
{		
	double d;
	gsl_vector *w=gsl_vector_alloc (G->size1);
	gsl_vector *Gw=gsl_vector_alloc (G->size1);
	gsl_vector_set_all (w, 1.0);
	
	gsl_blas_dgemv (CblasNoTrans, 1.0, G, w, 0.0, Gw);			
	gsl_blas_dsyr2 (CblasUpper, -1.0/(double)G->size1, Gw, w, G);
	gsl_blas_ddot (w, Gw, &d);		
	gsl_blas_dsyr (CblasUpper, d/((double)G->size1*(double)G->size1), w, G);
	
	for (size_t i=0; i<G->size1; ++i) {
		for (size_t j=0; j<i; ++j) {
			d=gsl_matrix_get (G, j, i);
			gsl_matrix_set (G, i, j, d);
		}
	}
	
	gsl_vector_free(w);
	gsl_vector_free(Gw);
	
	return;
}


//center the matrix G	
void CenterMatrix (gsl_matrix *G, gsl_vector *w)
{		
	double d, wtw;
	gsl_vector *Gw=gsl_vector_alloc (G->size1);
	
	gsl_blas_ddot (w, w, &wtw);	
	gsl_blas_dgemv (CblasNoTrans, 1.0, G, w, 0.0, Gw);			
	gsl_blas_dsyr2 (CblasUpper, -1.0/wtw, Gw, w, G);
	gsl_blas_ddot (w, Gw, &d);		
	gsl_blas_dsyr (CblasUpper, d/(wtw*wtw), w, G);
	
	for (size_t i=0; i<G->size1; ++i) {
		for (size_t j=0; j<i; ++j) {
			d=gsl_matrix_get (G, j, i);
			gsl_matrix_set (G, i, j, d);
		}
	}
	
	gsl_vector_free(Gw);
	
	return;
}


//scale the matrix G such that the mean diagonal = 1
void ScaleMatrix (gsl_matrix *G)
{		
	double d=0.0;
	
	for (size_t i=0; i<G->size1; ++i) {
		d+=gsl_matrix_get(G, i, i);
	}
	d/=(double)G->size1;
	
	gsl_matrix_scale (G, 1.0/d);
	
	return;
}


//center the vector y
double CenterVector (gsl_vector *y)
{		
	double d=0.0;
	
	for (size_t i=0; i<y->size; ++i) {
		d+=gsl_vector_get (y, i);
	}
	d/=(double)y->size;
	
	gsl_vector_add_constant (y, -1.0*d);
	
	return d;
}


//calculate UtX
void CalcUtX (const gsl_matrix *U, gsl_matrix *UtX) 
{
	gsl_vector *Utx_vec=gsl_vector_alloc (UtX->size1);
	for (size_t i=0; i<UtX->size2; ++i) {
		gsl_vector_view UtX_col=gsl_matrix_column (UtX, i);
		gsl_blas_dgemv (CblasTrans, 1.0, U, &UtX_col.vector, 0.0, Utx_vec);
		gsl_vector_memcpy (&UtX_col.vector, Utx_vec);
	}	
	gsl_vector_free (Utx_vec);
	return;
}


void CalcUtX (const gsl_matrix *U, const gsl_matrix *X, gsl_matrix *UtX) 
{
	for (size_t i=0; i<X->size2; ++i) {
		gsl_vector_const_view X_col=gsl_matrix_const_column (X, i);
		gsl_vector_view UtX_col=gsl_matrix_column (UtX, i);
		gsl_blas_dgemv (CblasTrans, 1.0, U, &X_col.vector, 0.0, &UtX_col.vector);
	}
	return;
}

void CalcUtX (const gsl_matrix *U, const gsl_vector *x, gsl_vector *Utx) 
{
	gsl_blas_dgemv (CblasTrans, 1.0, U, x, 0.0, Utx);
	return;
}

