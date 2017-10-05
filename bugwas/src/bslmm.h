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


#ifndef __BSLMM_H__                
#define __BSLMM_H__

#include <vector>
#include <map>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#ifdef FORCE_FLOAT
#include "param_float.h"
#else
#include "param.h"
#endif


using namespace std;






class BSLMM {

public:	
	// IO related parameters
	int a_mode;	
	size_t d_pace;
	
	string file_bfile;
	string file_geno;
	string file_out;
	
	// LMM related parameters
	double l_min;
	double l_max;
	size_t n_region;
	double pve_null;
	double pheno_mean;
	
	// BSLMM MCMC related parameters
	double h_min, h_max, h_scale;			//priors for h
	double rho_min, rho_max, rho_scale;		//priors for rho
	double logp_min, logp_max, logp_scale;		//priors for log(pi)
	size_t s_min, s_max;			//minimum and maximum number of gammas
	size_t w_step;					//number of warm up/burn in iterations
	size_t s_step;					//number of sampling iterations
	size_t r_pace;					//record pace
	size_t w_pace;					//write pace
	size_t n_accept;				//number of acceptance
	size_t n_mh;					//number of MH steps within each iteration
	double geo_mean;				//mean of the geometric distribution
	long int randseed;
	double trace_G;	
	
	HYPBSLMM cHyp_initial;

	// Summary statistics
	size_t ni_total, ns_total;	//number of total individuals and snps
	size_t ni_test, ns_test;	//number of individuals and snps used for analysis
	size_t n_cvt;				//number of covariates
	double time_UtZ;
	double time_Omega;		//time spent on optimization iterations
	
	vector<int> indicator_idv;				//indicator for individuals (phenotypes), 0 missing, 1 available for analysis
	vector<int> indicator_snp;				//sequence indicator for SNPs: 0 ignored because of (a) maf, (b) miss, (c) non-poly; 1 available for analysis
	
	vector<SNPINFO> snpInfo;		//record SNP information
	
	// Not included in PARAM
	gsl_rng *gsl_r; 
	gsl_ran_discrete_t *gsl_t;	
	map<size_t, size_t> mapRank2pos;	
	
	// Main Functions
	void CopyFromParam (PARAM &cPar);
	void CopyToParam (PARAM &cPar);
	
	void RidgeR(const gsl_matrix *U, const gsl_matrix *UtX, const gsl_vector *Uty, const gsl_vector *eval, const double lambda);
	
	void MCMC (const gsl_matrix *U, const gsl_matrix *UtX, const gsl_vector *Uty, const gsl_vector *K_eval, const gsl_vector *y);
	void WriteLog ();
	void WriteLR ();
	void WriteBV (const gsl_vector *bv);
	void WriteParam (vector<pair<double, double> > &beta_g, const gsl_vector *alpha, const size_t w);
	void WriteParam (const gsl_vector *alpha);
	void WriteResult (const int flag, const gsl_matrix *Result_hyp, const gsl_matrix *Result_gamma, const size_t w_col);
	
	//Subfunctions inside MCMC
	void CalcPgamma (double *p_gammar);
	
	double CalcPveLM (const gsl_matrix *UtXgamma, const gsl_vector *Uty, const double sigma_a2);
	void InitialMCMC (const gsl_matrix *UtX, const gsl_vector *Uty, vector<size_t> &rank_old, class HYPBSLMM &cHyp, vector<pair<size_t, double> > &pos_loglr);
	double CalcPosterior (const gsl_vector *Uty, const gsl_vector *K_eval, gsl_vector *Utu, gsl_vector *alpha_prime, class HYPBSLMM &cHyp);
	double CalcPosterior (const gsl_matrix *UtXgamma, const gsl_vector *Uty, const gsl_vector *K_eval, gsl_vector *UtXb, gsl_vector *Utu, gsl_vector *alpha_prime, gsl_vector *beta, class HYPBSLMM &cHyp);
	void CalcCC_PVEnZ (const gsl_matrix *U, const gsl_vector *Utu, gsl_vector *z_hat, class HYPBSLMM &cHyp);
	void CalcCC_PVEnZ (const gsl_matrix *U, const gsl_vector *UtXb, const gsl_vector *Utu, gsl_vector *z_hat, class HYPBSLMM &cHyp);
	double CalcREMLE (const gsl_matrix *Utw, const gsl_vector *Uty, const gsl_vector *K_eval);
	double CalcLR (const gsl_matrix *U, const gsl_matrix *UtX, const gsl_vector *Uty, const gsl_vector *K_eval, vector<pair<size_t, double> > &loglr_sort);		//calculate the maximum marginal likelihood ratio for each analyzed SNPs with gemma, use it to rank SNPs
	void SampleZ (const gsl_vector *y, const gsl_vector *z_hat, gsl_vector *z);
	double ProposeHnRho (const class HYPBSLMM &cHyp_old, class HYPBSLMM &cHyp_new, const size_t &repeat);
	double ProposePi (const class HYPBSLMM &cHyp_old, class HYPBSLMM &cHyp_new, const size_t &repeat);
	double ProposeGamma (const vector<size_t> &rank_old, vector<size_t> &rank_new, const double *p_gamma, const class HYPBSLMM &cHyp_old, class HYPBSLMM &cHyp_new, const size_t &repeat);
	void SetXgamma (gsl_matrix *Xgamma, const gsl_matrix *X, vector<size_t> &rank);
	
	//utility functions
//	double vec_sum (gsl_vector *v);
//	void vec_center (gsl_vector *v);
//	double calc_var (gsl_vector *v);
//	void calc_sigma (MCMC &cMcmc);
//	bool comp_lr (pair<size_t, double> a, pair<size_t, double> b);
};



#endif


