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

#ifndef __MATHFUNC_H__                
#define __MATHFUNC_H__

#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"


using namespace std;

double VectorVar (const gsl_vector *v);
void CenterMatrix (gsl_matrix *G);
void CenterMatrix (gsl_matrix *G, gsl_vector *w);
void ScaleMatrix (gsl_matrix *G);
double CenterVector (gsl_vector *y);
void CalcUtX (const gsl_matrix *U, gsl_matrix *UtX);
void CalcUtX (const gsl_matrix *U, const gsl_matrix *X, gsl_matrix *UtX);
void CalcUtX (const gsl_matrix *U, const gsl_vector *x, gsl_vector *Utx);



#endif
