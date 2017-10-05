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

#ifndef __LAPACK_H__                
#define __LAPACK_H__



using namespace std;


void lapack_float_cholesky_decomp (gsl_matrix_float *A);
void lapack_cholesky_decomp (gsl_matrix *A);
void lapack_float_cholesky_solve (gsl_matrix_float *A, const gsl_vector_float *b, gsl_vector_float *x);
void lapack_cholesky_solve (gsl_matrix *A, const gsl_vector *b, gsl_vector *x);
void lapack_sgemm (char *TransA, char *TransB, float alpha, const gsl_matrix_float *A, const gsl_matrix_float *B, float beta, gsl_matrix_float *C);
void lapack_dgemm (char *TransA, char *TransB, double alpha, const gsl_matrix *A, const gsl_matrix *B, double beta, gsl_matrix *C);
void lapack_float_eigen_symmv (gsl_matrix_float *A, gsl_vector_float *eval, gsl_matrix_float *evec);
void lapack_eigen_symmv (gsl_matrix *A, gsl_vector *eval, gsl_matrix *evec);

double EigenDecomp (gsl_matrix *G, gsl_matrix *U, gsl_vector *eval);
double EigenDecomp (gsl_matrix_float *G, gsl_matrix_float *U, gsl_vector_float *eval);

double CholeskySolve(gsl_matrix *Omega, gsl_vector *Xty, gsl_vector *OiXty);
double CholeskySolve(gsl_matrix_float *Omega, gsl_vector_float *Xty, gsl_vector_float *OiXty);

#endif



