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

#ifndef __PRDT_H__                
#define __PRDT_H__


#include <vector>
#include <map>
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"

#ifdef FORCE_FLOAT
#include "param_float.h"
#else
#include "param.h"
#endif

using namespace std;

class PRDT {
	
public:
	// IO related parameters
	size_t a_mode;
	size_t d_pace;
	
	string file_bfile;
	string file_geno;
	string file_out;
	
	vector<int> indicator_idv;
	vector<SNPINFO> snpInfo;
	map<string, double> mapRS2est;
	
	size_t ns_total;
	size_t ns_test;
	
	double time_eigen;
	
	// Main functions
	void CopyFromParam (PARAM &cPar);
	void CopyToParam (PARAM &cPar);
	void WriteFiles (gsl_vector *y_prdt);
	void AddBV (gsl_matrix *G, const gsl_vector *u_hat, gsl_vector *y_prdt);
	void AnalyzeBimbam (gsl_vector *y_prdt);
	void AnalyzePlink (gsl_vector *y_prdt);
};


#endif







