#include <iostream>
#include <fstream>
#include <sstream>

#include <iomanip>
#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h> 
#include <ctime>
#include <cstring>
#include <algorithm>

#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_eigen.h"
#include "gsl/gsl_randist.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_roots.h"




//#ifdef WITH_LAPACK
#include "lapack.h"
//#endif

#ifdef FORCE_FLOAT
#include "param_float.h"
#include "bslmm_float.h"
#include "lmm_float.h"  //for class FUNC_PARAM and MatrixCalcLR
#include "mathfunc_float.h"  //for function CenterVector
#else
#include "param.h"
#include "bslmm.h"
#include "lmm.h"
#include "mathfunc.h"
#endif

using namespace std;




void BSLMM::CopyFromParam (PARAM &cPar) 
{
	a_mode=cPar.a_mode;
	d_pace=cPar.d_pace;
	
	file_bfile=cPar.file_bfile;
	file_geno=cPar.file_geno;
	file_out=cPar.file_out;
	
	l_min=cPar.h_min;	
	l_max=cPar.h_max;  
	n_region=cPar.n_region;	
	pve_null=cPar.pve_null;
	pheno_mean=cPar.pheno_mean;
	
	time_UtZ=0.0;
	time_Omega=0.0;
	n_accept=0;
	
	h_min=cPar.h_min;	
	h_max=cPar.h_max;  
	h_scale=cPar.h_scale;
	rho_min=cPar.rho_min;	
	rho_max=cPar.rho_max;  
	rho_scale=cPar.rho_scale;
	logp_min=cPar.logp_min;	
	logp_max=cPar.logp_max;  
	logp_scale=cPar.logp_scale;
	
	s_min=cPar.s_min;
	s_max=cPar.s_max;
	w_step=cPar.w_step;
	s_step=cPar.s_step;
	r_pace=cPar.r_pace;
	w_pace=cPar.w_pace;
	n_mh=cPar.n_mh;
	geo_mean=cPar.geo_mean;
	randseed=cPar.randseed;
	trace_G=cPar.trace_G;
	
	ni_total=cPar.ni_total;
	ns_total=cPar.ns_total;
	ni_test=cPar.ni_test;
	ns_test=cPar.ns_test;
	n_cvt=cPar.n_cvt;
	
	indicator_idv=cPar.indicator_idv;
	indicator_snp=cPar.indicator_snp;
	snpInfo=cPar.snpInfo;
	
	return;
}


void BSLMM::CopyToParam (PARAM &cPar) 
{
	cPar.time_UtZ=time_UtZ;
	cPar.time_Omega=time_Omega;
	cPar.cHyp_initial=cHyp_initial;
	cPar.n_accept=n_accept;
	cPar.pheno_mean=pheno_mean;
	cPar.randseed=randseed;
	
	return;
}



void BSLMM::WriteBV (const gsl_vector *bv) 
{
	string file_str;
	file_str="./output/"+file_out;
	file_str+=".bv.txt";

	ofstream outfile (file_str.c_str(), ofstream::out);
	if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
	
	size_t t=0;
	for (size_t i=0; i<ni_total; ++i) {
		if (indicator_idv[i]==0) {
			outfile<<"NA"<<endl;
		}		
		else {
			outfile<<scientific<<setprecision(6)<<gsl_vector_get(bv, t)<<endl;
			t++;
		}
	}		
	
	outfile.clear();	
	outfile.close();	
	return;
}




void BSLMM::WriteParam (vector<pair<double, double> > &beta_g, const gsl_vector *alpha, const size_t w) 
{
	string file_str;
	file_str="./output/"+file_out;
	file_str+=".param.txt";

	ofstream outfile (file_str.c_str(), ofstream::out);
	if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
	
	outfile<<"chr"<<"\t"<<"rs"<<"\t"
			<<"ps"<<"\t"<<"n_miss"<<"\t"<<"alpha"<<"\t"
			<<"beta"<<"\t"<<"gamma"<<endl;
	
	size_t t=0;
	for (size_t i=0; i<ns_total; ++i) {
		if (indicator_snp[i]==0) {continue;}		
		
		outfile<<snpInfo[i].chr<<"\t"<<snpInfo[i].rs_number<<"\t"
		<<snpInfo[i].base_position<<"\t"<<snpInfo[i].n_miss<<"\t";	
				
		outfile<<scientific<<setprecision(6)<<gsl_vector_get(alpha, t)<<"\t";
		if (beta_g[t].second!=0) {
			outfile<<beta_g[t].first/beta_g[t].second<<"\t"<<beta_g[t].second/(double)w<<endl;
		}
		else {
			outfile<<0.0<<"\t"<<0.0<<endl;
		}
		t++;
	}		
	
	outfile.clear();	
	outfile.close();	
	return;
}


void BSLMM::WriteParam (const gsl_vector *alpha) 
{
	string file_str;
	file_str="./output/"+file_out;
	file_str+=".param.txt";

	ofstream outfile (file_str.c_str(), ofstream::out);
	if (!outfile) {cout<<"error writing file: "<<file_str.c_str()<<endl; return;}
	
	outfile<<"chr"<<"\t"<<"rs"<<"\t"
			<<"ps"<<"\t"<<"n_miss"<<"\t"<<"alpha"<<"\t"
			<<"beta"<<"\t"<<"gamma"<<endl;
	
	size_t t=0;
	for (size_t i=0; i<ns_total; ++i) {
		if (indicator_snp[i]==0) {continue;}		

		outfile<<snpInfo[i].chr<<"\t"<<snpInfo[i].rs_number<<"\t"
				<<snpInfo[i].base_position<<"\t"<<snpInfo[i].n_miss<<"\t";				
		outfile<<scientific<<setprecision(6)<<gsl_vector_get(alpha, t)<<"\t";
		outfile<<0.0<<"\t"<<0.0<<endl;
		t++;
	}		
	
	outfile.clear();	
	outfile.close();	
	return;
}


void BSLMM::WriteResult (const int flag, const gsl_matrix *Result_hyp, const gsl_matrix *Result_gamma, const size_t w_col) 
{
	string file_gamma, file_hyp;
	file_gamma="./output/"+file_out;
	file_gamma+=".gamma.txt";
	file_hyp="./output/"+file_out;
	file_hyp+=".hyp.txt";

	ofstream outfile_gamma, outfile_hyp;
		
	if (flag==0) {
		outfile_gamma.open (file_gamma.c_str(), ofstream::out);
		outfile_hyp.open (file_hyp.c_str(), ofstream::out);
		if (!outfile_gamma) {cout<<"error writing file: "<<file_gamma<<endl; return;}
		if (!outfile_hyp) {cout<<"error writing file: "<<file_hyp<<endl; return;}
		
		outfile_hyp<<"h \t pve \t rho \t pge \t pi \t n_gamma"<<endl;
		
		for (size_t i=0; i<s_max; ++i) {
			outfile_gamma<<"s"<<i<<"\t";
		}
		outfile_gamma<<endl;
	}
	else {
		outfile_gamma.open (file_gamma.c_str(), ofstream::app);
		outfile_hyp.open (file_hyp.c_str(), ofstream::app);
		if (!outfile_gamma) {cout<<"error writing file: "<<file_gamma<<endl; return;}
		if (!outfile_hyp) {cout<<"error writing file: "<<file_hyp<<endl; return;}
		
		size_t w;
		if (w_col==0) {w=w_pace;}
		else {w=w_col;}
		
		for (size_t i=0; i<w; ++i) {
			outfile_hyp<<scientific;
			for (size_t j=0; j<4; ++j) {
				outfile_hyp<<setprecision(6)<<gsl_matrix_get (Result_hyp, i, j)<<"\t";
			}
			outfile_hyp<<setprecision(6)<<exp(gsl_matrix_get (Result_hyp, i, 4))<<"\t";
			outfile_hyp<<(int)gsl_matrix_get (Result_hyp, i, 5)<<"\t";
			outfile_hyp<<endl;
		}
		
		for (size_t i=0; i<w; ++i) {
			for (size_t j=0; j<s_max; ++j) {
				outfile_gamma<<(int)gsl_matrix_get (Result_gamma, i, j)<<"\t";
			}
			outfile_gamma<<endl;
		}
		
	}
	
	outfile_hyp.close();
	outfile_hyp.clear();
	outfile_gamma.close();
	outfile_gamma.clear();	
	return;
}



void BSLMM::CalcPgamma (double *p_gamma)
{
	double p, s=0.0;
	for (size_t i=0; i<ns_test; ++i) {
		p=0.7*gsl_ran_geometric_pdf (i+1, 1.0/geo_mean)+0.3/(double)ns_test;
		p_gamma[i]=p;
		s+=p;
	}
	for (size_t i=0; i<ns_test; ++i) {
		p=p_gamma[i];
		p_gamma[i]=p/s;
	}
	return;
}



void BSLMM::SetXgamma (gsl_matrix *Xgamma, const gsl_matrix *X, vector<size_t> &rank)
{
	size_t pos;
	for (size_t i=0; i<rank.size(); ++i) {
		pos=mapRank2pos[rank[i]];
		gsl_vector_view Xgamma_col=gsl_matrix_column (Xgamma, i);
		gsl_vector_const_view X_col=gsl_matrix_const_column (X, pos);
		gsl_vector_memcpy (&Xgamma_col.vector, &X_col.vector);
	}
	
	return;
}



double BSLMM::CalcPveLM (const gsl_matrix *UtXgamma, const gsl_vector *Uty, const double sigma_a2) 
{
	double pve, var_y;	
	
	gsl_matrix *Omega=gsl_matrix_alloc (UtXgamma->size2, UtXgamma->size2);
	gsl_vector *Xty=gsl_vector_alloc (UtXgamma->size2);
	gsl_vector *OiXty=gsl_vector_alloc (UtXgamma->size2);

	gsl_matrix_set_identity (Omega);
	gsl_matrix_scale (Omega, 1.0/sigma_a2); 

#ifdef WITH_LAPACK
	lapack_dgemm ((char *)"T", (char *)"N", 1.0, UtXgamma, UtXgamma, 1.0, Omega);
#else
	gsl_blas_dgemm (CblasTrans, CblasNoTrans, 1.0, UtXgamma, UtXgamma, 1.0, Omega);	
#endif
	gsl_blas_dgemv (CblasTrans, 1.0, UtXgamma, Uty, 0.0, Xty);

	CholeskySolve(Omega, Xty, OiXty);
	
	gsl_blas_ddot (Xty, OiXty, &pve);
	gsl_blas_ddot (Uty, Uty, &var_y);
	
	pve/=var_y;
	
	gsl_matrix_free (Omega);
	gsl_vector_free (Xty);
	gsl_vector_free (OiXty);

	return pve;
}


void BSLMM::InitialMCMC (const gsl_matrix *UtX, const gsl_vector *Uty, vector<size_t> &rank, class HYPBSLMM &cHyp, vector<pair<size_t, double> > &pos_loglr)
{
	double q_genome=gsl_cdf_chisq_Qinv(0.05/(double)ns_test, 1);
	
	cHyp.n_gamma=0;
	for (size_t i=0; i<pos_loglr.size(); ++i) {
		if (2.0*pos_loglr[i].second>q_genome) {cHyp.n_gamma++;}
	}
	if (cHyp.n_gamma<10) {cHyp.n_gamma=10;}
	
	if (cHyp.n_gamma>s_max) {cHyp.n_gamma=s_max;}
	if (cHyp.n_gamma<s_min) {cHyp.n_gamma=s_min;}	
	
	rank.clear();
	for (size_t i=0; i<cHyp.n_gamma; ++i) {
		rank.push_back(i);
	}
	
	cHyp.logp=log((double)cHyp.n_gamma/(double)ns_test);
	cHyp.h=pve_null; 
	
	if (cHyp.logp==0) {cHyp.logp=-0.000001;}
	
	gsl_matrix *UtXgamma=gsl_matrix_alloc (ni_test, cHyp.n_gamma);
	SetXgamma (UtXgamma, UtX, rank);
	double sigma_a2=cHyp.h*1.0/(trace_G*(1-cHyp.h)*exp(cHyp.logp)*(double)ns_test);
	if (sigma_a2==0) {sigma_a2=0.025;}	
	cHyp.rho=CalcPveLM (UtXgamma, Uty, sigma_a2)/cHyp.h;
	gsl_matrix_free (UtXgamma);
	
	if (cHyp.rho>1.0) {cHyp.rho=1.0;}
	
	if (cHyp.h<h_min) {cHyp.h=h_min;}
	if (cHyp.h>h_max) {cHyp.h=h_max;}
	if (cHyp.rho<rho_min) {cHyp.rho=rho_min;}
	if (cHyp.rho>rho_max) {cHyp.rho=rho_max;}
	if (cHyp.logp<logp_min) {cHyp.logp=logp_min;}
	if (cHyp.logp>logp_max) {cHyp.logp=logp_max;}
	
	
//	if (fix_sigma>=0) {
//		fix_sigma=cHyp.h;
//		rho_max=1-cHyp.h;
//		cHyp.rho=rho_max/2.0;
//	}
	
	//Initial for grid sampling:
//	cHyp.h=0.225;
//	cHyp.rho=1.0;
//	cHyp.logp=-4.835429;
	
	cout<<"initial value of h = "<<cHyp.h<<endl;
	cout<<"initial value of rho = "<<cHyp.rho<<endl;
	cout<<"initial value of pi = "<<exp(cHyp.logp)<<endl;
	cout<<"initial value of |gamma| = "<<cHyp.n_gamma<<endl;
	
	return;
}



double BSLMM::CalcPosterior (const gsl_vector *Uty, const gsl_vector *K_eval, gsl_vector *Utu, gsl_vector *alpha_prime, class HYPBSLMM &cHyp)
{
	double sigma_b2=cHyp.h*(1.0-cHyp.rho)/(trace_G*(1-cHyp.h));
	
	gsl_vector *Utu_rand=gsl_vector_alloc (Uty->size);	
	gsl_vector *weight_Hi=gsl_vector_alloc (Uty->size);
	
	double logpost=0.0;
	double d, ds, uy, Hi_yy=0, logdet_H=0.0;
	for (size_t i=0; i<ni_test; ++i) {
		d=gsl_vector_get (K_eval, i)*sigma_b2;
		ds=d/(d+1.0);
		d=1.0/(d+1.0);		
		gsl_vector_set (weight_Hi, i, d);
		
		logdet_H-=log(d);
		uy=gsl_vector_get (Uty, i);
		Hi_yy+=d*uy*uy;
		
		gsl_vector_set (Utu_rand, i, gsl_ran_gaussian(gsl_r, 1)*sqrt(ds));
	}
	
	//sample tau
	double tau=1.0;
	if (a_mode==11) {tau = gsl_ran_gamma (gsl_r, (double)ni_test/2.0,  2.0/Hi_yy); }
	
	//sample alpha
	gsl_vector_memcpy (alpha_prime, Uty);
	gsl_vector_mul (alpha_prime, weight_Hi);
	gsl_vector_scale (alpha_prime, sigma_b2);
	
	//sample u
	gsl_vector_memcpy (Utu, alpha_prime);
	gsl_vector_mul (Utu, K_eval);
	if (a_mode==11) {gsl_vector_scale (Utu_rand, sqrt(1.0/tau));}
	gsl_vector_add (Utu, Utu_rand);	
	
	//for quantitative traits, calculate pve and ppe
	if (a_mode==11) {
		gsl_blas_ddot (Utu, Utu, &d);
		cHyp.pve=d/(double)ni_test;	
		cHyp.pve/=cHyp.pve+1.0/tau;
		cHyp.pge=0.0;	
	}

	//calculate likelihood
	logpost=-0.5*logdet_H;
	if (a_mode==11) {logpost-=0.5*(double)ni_test*log(Hi_yy);}
	else {logpost-=0.5*Hi_yy;}
	
	logpost+=((double)cHyp.n_gamma-1.0)*cHyp.logp+((double)ns_test-(double)cHyp.n_gamma)*log(1-exp(cHyp.logp));
	
	gsl_vector_free (Utu_rand);
	gsl_vector_free (weight_Hi);
	
	return logpost;
}


double BSLMM::CalcPosterior (const gsl_matrix *UtXgamma, const gsl_vector *Uty, const gsl_vector *K_eval, gsl_vector *UtXb, gsl_vector *Utu, gsl_vector *alpha_prime, gsl_vector *beta, class HYPBSLMM &cHyp)
{
	clock_t time_start;	
	
	double sigma_a2=cHyp.h*cHyp.rho/(trace_G*(1-cHyp.h)*exp(cHyp.logp)*(double)ns_test);
	double sigma_b2=cHyp.h*(1.0-cHyp.rho)/(trace_G*(1-cHyp.h));
	
	double logpost=0.0;
	double d, ds, uy, P_yy=0, logdet_O=0.0, logdet_H=0.0;
	
	gsl_matrix *UtXgamma_eval=gsl_matrix_alloc (UtXgamma->size1, UtXgamma->size2);	
	gsl_matrix *Omega=gsl_matrix_alloc (UtXgamma->size2, UtXgamma->size2);
	gsl_vector *XtHiy=gsl_vector_alloc (UtXgamma->size2);
	gsl_vector *beta_hat=gsl_vector_alloc (UtXgamma->size2);
	gsl_vector *Utu_rand=gsl_vector_alloc (UtXgamma->size1);	
	gsl_vector *weight_Hi=gsl_vector_alloc (UtXgamma->size1);
	
	gsl_matrix_memcpy (UtXgamma_eval, UtXgamma);
	
	logdet_H=0.0; P_yy=0.0;
	for (size_t i=0; i<ni_test; ++i) {
		gsl_vector_view UtXgamma_row=gsl_matrix_row (UtXgamma_eval, i);
		d=gsl_vector_get (K_eval, i)*sigma_b2;
		ds=d/(d+1.0);
		d=1.0/(d+1.0);
		gsl_vector_set (weight_Hi, i, d);
		
		logdet_H-=log(d);
		uy=gsl_vector_get (Uty, i);
		P_yy+=d*uy*uy;
		gsl_vector_scale (&UtXgamma_row.vector, d);
		
		gsl_vector_set (Utu_rand, i, gsl_ran_gaussian(gsl_r, 1)*sqrt(ds));
	}
	
	//calculate Omega
	gsl_matrix_set_identity (Omega);
	
	time_start=clock();
#ifdef WITH_LAPACK
	lapack_dgemm ((char *)"T", (char *)"N", sigma_a2, UtXgamma_eval, UtXgamma, 1.0, Omega);
#else
	gsl_blas_dgemm (CblasTrans, CblasNoTrans, sigma_a2, UtXgamma_eval, UtXgamma, 1.0, Omega);
#endif	
	time_Omega+=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);
	
	
	//calculate beta_hat
	gsl_blas_dgemv (CblasTrans, 1.0, UtXgamma_eval, Uty, 0.0, XtHiy);	

	logdet_O=CholeskySolve(Omega, XtHiy, beta_hat);
	
	gsl_vector_scale (beta_hat, sigma_a2);

	gsl_blas_ddot (XtHiy, beta_hat, &d);
	P_yy-=d;
	
	//sample tau
	double tau=1.0;
	if (a_mode==11) {tau =gsl_ran_gamma (gsl_r, (double)ni_test/2.0,  2.0/P_yy); }

	//sample beta
	for (size_t i=0; i<beta->size; i++)
	{
		d=gsl_ran_gaussian(gsl_r, 1); 
		gsl_vector_set(beta, i, d); 
	}
	gsl_blas_dtrsv(CblasUpper, CblasNoTrans, CblasNonUnit, Omega, beta); 
	
	
	//it compuates inv(L^T(Omega)) %*% beta;  
	gsl_vector_scale(beta, sqrt(sigma_a2/tau));
	gsl_vector_add(beta, beta_hat); 
	gsl_blas_dgemv (CblasNoTrans, 1.0, UtXgamma, beta, 0.0, UtXb);
	
	//sample alpha
	gsl_vector_memcpy (alpha_prime, Uty);
	gsl_vector_sub (alpha_prime, UtXb);
	gsl_vector_mul (alpha_prime, weight_Hi);
	gsl_vector_scale (alpha_prime, sigma_b2);
	
	//sample u
	gsl_vector_memcpy (Utu, alpha_prime);
	gsl_vector_mul (Utu, K_eval);
	
	if (a_mode==11) {gsl_vector_scale (Utu_rand, sqrt(1.0/tau));}
	gsl_vector_add (Utu, Utu_rand);	
	
	
	//for quantitative traits, calculate pve and pge
	if (a_mode==11) {
		gsl_blas_ddot (UtXb, UtXb, &d);
		cHyp.pge=d/(double)ni_test;
	
		gsl_blas_ddot (Utu, Utu, &d);
		cHyp.pve=cHyp.pge+d/(double)ni_test;
		
		if (cHyp.pve==0) {cHyp.pge=0.0;}
		else {cHyp.pge/=cHyp.pve;}
		cHyp.pve/=cHyp.pve+1.0/tau;	
	}	
		
	gsl_matrix_free (UtXgamma_eval);
	gsl_matrix_free (Omega);
	gsl_vector_free (XtHiy);
	gsl_vector_free (beta_hat);
	gsl_vector_free (Utu_rand);	
	gsl_vector_free (weight_Hi);
	
	logpost=-0.5*logdet_H-0.5*logdet_O;
	if (a_mode==11) {logpost-=0.5*(double)ni_test*log(P_yy);}
	else {logpost-=0.5*P_yy;}
//	else {logpost+=-0.5*P_yy*tau+0.5*(double)ni_test*log(tau);}
	logpost+=((double)cHyp.n_gamma-1.0)*cHyp.logp+((double)ns_test-(double)cHyp.n_gamma)*log(1.0-exp(cHyp.logp));
	
	return logpost;
}



//calculate pve and pge, and calculate z_hat for case-control data	
void BSLMM::CalcCC_PVEnZ (const gsl_matrix *U, const gsl_vector *Utu, gsl_vector *z_hat, class HYPBSLMM &cHyp) 
{
	double d;
	
	gsl_blas_ddot (Utu, Utu, &d);
	cHyp.pve=d/(double)ni_test;	
		
	gsl_blas_dgemv (CblasNoTrans, 1.0, U, Utu, 0.0, z_hat);
		
	cHyp.pve/=cHyp.pve+1.0;
	cHyp.pge=0.0;	
	
	return;
}


//calculate pve and pge, and calculate z_hat for case-control data	
void BSLMM::CalcCC_PVEnZ (const gsl_matrix *U, const gsl_vector *UtXb, const gsl_vector *Utu, gsl_vector *z_hat, class HYPBSLMM &cHyp) 
{
	double d;
	gsl_vector *UtXbU=gsl_vector_alloc (Utu->size);
	
	gsl_blas_ddot (UtXb, UtXb, &d);
	cHyp.pge=d/(double)ni_test;
	
	gsl_blas_ddot (Utu, Utu, &d);
	cHyp.pve=cHyp.pge+d/(double)ni_test;
	
	gsl_vector_memcpy (UtXbU, Utu);
	gsl_vector_add (UtXbU, UtXb);
	gsl_blas_dgemv (CblasNoTrans, 1.0, U, UtXbU, 0.0, z_hat);	
	
	if (cHyp.pve==0) {cHyp.pge=0.0;}
	else {cHyp.pge/=cHyp.pve;}
	
	cHyp.pve/=cHyp.pve+1.0;
	
	gsl_vector_free(UtXbU);
	return;
}




void BSLMM::SampleZ (const gsl_vector *y, const gsl_vector *z_hat, gsl_vector *z)
{	
	double d1, d2, z_rand=0.0;
	for (size_t i=0; i<z->size; ++i) {
		d1=gsl_vector_get (y, i);
		d2=gsl_vector_get (z_hat, i);
		//0.1 is chosen arbitrary, any small number will be fine
		if (d1<0.1) {
			//control, right truncated
			do {				
				z_rand=d2+gsl_ran_gaussian(gsl_r, 1.0);
			} while (z_rand>0.0);
		}
		else {
			do {
				z_rand=d2+gsl_ran_gaussian(gsl_r, 1.0);
			} while (z_rand<0.0);
		}
		
		gsl_vector_set (z, i, z_rand);
	}

	return;
}





double BSLMM::ProposeHnRho (const class HYPBSLMM &cHyp_old, class HYPBSLMM &cHyp_new, const size_t &repeat)
{
	
	double h=cHyp_old.h, rho=cHyp_old.rho;
	
	double d_h=(h_max-h_min)*h_scale, d_rho=(rho_max-rho_min)*rho_scale;
	
	for (size_t i=0; i<repeat; ++i) {
		h=h+(gsl_rng_uniform(gsl_r)-0.5)*d_h;
		if (h<h_min) {h=2*h_min-h;}
		if (h>h_max) {h=2*h_max-h;}
		
		rho=rho+(gsl_rng_uniform(gsl_r)-0.5)*d_rho;
		if (rho<rho_min) {rho=2*rho_min-rho;}
		if (rho>rho_max) {rho=2*rho_max-rho;}
	}
	/*
	//Grid Sampling
	for (size_t i=0; i<repeat; ++i) {
		if (gsl_rng_uniform(gsl_r)<0.66) {continue;}
		h=h+(gsl_rng_uniform_int(gsl_r, 2)-0.5)*0.1;
		if (h<h_min) {h=h_max;}
		if (h>h_max) {h=h_min;}
	}
	
	for (size_t i=0; i<repeat; ++i) {
		if (gsl_rng_uniform(gsl_r)<0.66) {continue;}
		rho=rho+(gsl_rng_uniform_int(gsl_r, 2)-0.5)*0.1;
		if (rho<rho_min) {rho=rho_max;}
		if (rho>rho_max) {rho=rho_min;}
	}
	*/
	cHyp_new.h=h;
	cHyp_new.rho=rho;
	return 0.0;
}


double BSLMM::ProposePi (const class HYPBSLMM &cHyp_old, class HYPBSLMM &cHyp_new, const size_t &repeat)
{
	double logp_old=cHyp_old.logp, logp_new=cHyp_old.logp;
	double log_ratio=0.0;
	
	double d_logp=min(0.1, (logp_max-logp_min)*logp_scale);
	
	for (size_t i=0; i<repeat; ++i) {
		logp_new=logp_old+(gsl_rng_uniform(gsl_r)-0.5)*d_logp;
		if (logp_new<logp_min) {logp_new=2*logp_min-logp_new;}
		if (logp_new>logp_max) {logp_new=2*logp_max-logp_new;}		
		
		log_ratio+=logp_new-logp_old;
		logp_old=logp_new;
	}
	/*
	//Grid Sampling
	for (size_t i=0; i<repeat; ++i) {
		if (gsl_rng_uniform(gsl_r)<0.66) {continue;}
		logp_new=logp_old+(gsl_rng_uniform_int(gsl_r, 2)-0.5)*0.5*log(10.0);
		if (logp_new<logp_min) {logp_new=logp_max;}
		if (logp_new>logp_max) {logp_new=logp_min;}	
		
		log_ratio+=logp_new-logp_old;
		logp_old=logp_new;
	}
	*/
	cHyp_new.logp=logp_new;
	
	return log_ratio;
}

bool comp_vec (size_t a, size_t b)
{
	return (a < b); 
}


double BSLMM::ProposeGamma (const vector<size_t> &rank_old, vector<size_t> &rank_new, const double *p_gamma, const class HYPBSLMM &cHyp_old, class HYPBSLMM &cHyp_new, const size_t &repeat)
{
	map<size_t, int> mapRank2in;
	size_t r;
	double unif, logp=0.0;
	int flag_gamma;
	size_t r_add, r_remove, col_id;
	
	rank_new.clear();
	if (cHyp_old.n_gamma!=rank_old.size()) {cout<<"size wrong"<<endl;}
	
	if (cHyp_old.n_gamma!=0) {
		for (size_t i=0; i<rank_old.size(); ++i) {
			r=rank_old[i];
			rank_new.push_back(r);
			mapRank2in[r]=1;
		}
	}
	cHyp_new.n_gamma=cHyp_old.n_gamma;	
	
	for (size_t i=0; i<repeat; ++i) {
		unif=gsl_rng_uniform(gsl_r); 
	
		if (unif < 0.40 && cHyp_new.n_gamma<s_max) {flag_gamma=1;}
		else if (unif>=0.40 && unif < 0.80 && cHyp_new.n_gamma>s_min) {flag_gamma=2;}
		else if (unif>=0.80 && cHyp_new.n_gamma>0 && cHyp_new.n_gamma<ns_test) {flag_gamma=3;}
		else {flag_gamma=4;}
	
		if(flag_gamma==1)  {//add a snp; 
			do {
				r_add=gsl_ran_discrete (gsl_r, gsl_t);
			} while (mapRank2in.count(r_add)!=0); 
		
			double prob_total=1.0;
			for (size_t i=0; i<cHyp_new.n_gamma; ++i) {
				r=rank_new[i];
				prob_total-=p_gamma[r];
			}

			mapRank2in[r_add]=1;
			rank_new.push_back(r_add);
			cHyp_new.n_gamma++;
			logp+=-log(p_gamma[r_add]/prob_total)-log((double)cHyp_new.n_gamma);
		}
		else if (flag_gamma==2) {//delete a snp;
			col_id=gsl_rng_uniform_int(gsl_r, cHyp_new.n_gamma);		
			r_remove=rank_new[col_id];
		
			double prob_total=1.0;
			for (size_t i=0; i<cHyp_new.n_gamma; ++i) {
				r=rank_new[i];
				prob_total-=p_gamma[r];
			}
			prob_total+=p_gamma[r_remove];
		
			mapRank2in.erase(r_remove);
			rank_new.erase(rank_new.begin()+col_id);
			logp+=log(p_gamma[r_remove]/prob_total)+log((double)cHyp_new.n_gamma);
			cHyp_new.n_gamma--;
		}
		else if (flag_gamma==3) {//switch a snp;
			col_id=gsl_rng_uniform_int(gsl_r, cHyp_new.n_gamma);		
			r_remove=rank_new[col_id];
		//careful with the proposal
			do {
				r_add=gsl_ran_discrete (gsl_r, gsl_t);
			} while (mapRank2in.count(r_add)!=0); 
			
			double prob_total=1.0;
			for (size_t i=0; i<cHyp_new.n_gamma; ++i) {
				r=rank_new[i];
				prob_total-=p_gamma[r];
			}
			
			logp+=log(p_gamma[r_remove]/(prob_total+p_gamma[r_remove]-p_gamma[r_add]) );
			logp-=log(p_gamma[r_add]/prob_total);
			
			mapRank2in.erase(r_remove);
			mapRank2in[r_add]=1;
			rank_new.erase(rank_new.begin()+col_id);
			rank_new.push_back(r_add);
		}
		else {logp+=0;}//do not change
	}
	
	stable_sort (rank_new.begin(), rank_new.end(), comp_vec);

	mapRank2in.clear();
	return logp;
}






bool comp_lr (pair<size_t, double> a, pair<size_t, double> b)
{
	return (a.second > b.second); 
}







//if a_mode==13 then Uty==y
void BSLMM::MCMC (const gsl_matrix *U, const gsl_matrix *UtX, const gsl_vector *Uty, const gsl_vector *K_eval, const gsl_vector *y) {
	clock_t time_start;	

	class HYPBSLMM cHyp_old, cHyp_new;
	
	gsl_matrix *Result_hyp=gsl_matrix_alloc (w_pace, 6);
	gsl_matrix *Result_gamma=gsl_matrix_alloc (w_pace, s_max);	
	
	gsl_vector *alpha_prime=gsl_vector_alloc (ni_test);		
	gsl_vector *alpha_new=gsl_vector_alloc (ni_test);
	gsl_vector *alpha_old=gsl_vector_alloc (ni_test);	
	gsl_vector *Utu=gsl_vector_alloc (ni_test);
	gsl_vector *Utu_new=gsl_vector_alloc (ni_test);
	gsl_vector *Utu_old=gsl_vector_alloc (ni_test);
	
	gsl_vector *UtXb_new=gsl_vector_alloc (ni_test);
	gsl_vector *UtXb_old=gsl_vector_alloc (ni_test);
	
	gsl_vector *z_hat=gsl_vector_alloc (ni_test);
	gsl_vector *z=gsl_vector_alloc (ni_test);
	gsl_vector *Utz=gsl_vector_alloc (ni_test);	

	gsl_vector_memcpy (Utz, Uty);			
	
	double logPost_new, logPost_old;
	double logMHratio;
	double mean_z=0.0;	
	
	gsl_matrix_set_zero (Result_gamma);
	gsl_vector_set_zero (Utu);
	gsl_vector_set_zero (alpha_prime);
	if (a_mode==13) {
		pheno_mean=0.0;
	}
	
	vector<pair<double, double> > beta_g;
	for (size_t i=0; i<ns_test; i++) {
		beta_g.push_back(make_pair(0.0, 0.0));
	}
	
	vector<size_t> rank_new, rank_old;
	vector<double> beta_new, beta_old;	

	vector<pair<size_t, double> > pos_loglr;	
	
	MatrixCalcLR (U, UtX, Utz, K_eval, l_min, l_max, n_region, pos_loglr);

	stable_sort (pos_loglr.begin(), pos_loglr.end(), comp_lr);
	for (size_t i=0; i<ns_test; ++i) {
		mapRank2pos[i]=pos_loglr[i].first;
	}
	
	//calculate proposal distribution for gamma (unnormalized), and set up gsl_r and gsl_t			
	gsl_rng_env_setup();                
	const gsl_rng_type * gslType;                                               
	gslType = gsl_rng_default; 
	if (randseed<0)
	{
		time_t rawtime;
		time (&rawtime);
		tm * ptm = gmtime (&rawtime);
		
		randseed = (unsigned) (ptm->tm_hour%24*3600+ptm->tm_min*60+ptm->tm_sec);
	}
	gsl_r = gsl_rng_alloc(gslType); 
	gsl_rng_set(gsl_r, randseed);
	
	double *p_gamma = new double[ns_test]; 
	CalcPgamma (p_gamma);
	
	gsl_t=gsl_ran_discrete_preproc (ns_test, p_gamma);
	
	//initial parameters
	InitialMCMC (UtX, Utz, rank_old, cHyp_old, pos_loglr);
//	if (fix_sigma>=0) {
//		rho_max=1-fix_sigma;
//		cHyp_old.h=fix_sigma/(1-cHyp_old.rho);
//	}
	
	cHyp_initial=cHyp_old;
	
	if (cHyp_old.n_gamma==0 || cHyp_old.rho==0) {
		logPost_old=CalcPosterior(Utz, K_eval, Utu_old, alpha_old, cHyp_old);

		beta_old.clear();
		for (size_t i=0; i<cHyp_old.n_gamma; ++i) {
		  beta_old.push_back(0);
		}	
	}
	else {
		gsl_matrix *UtXgamma=gsl_matrix_alloc (ni_test, cHyp_old.n_gamma);
		gsl_vector *beta=gsl_vector_alloc (cHyp_old.n_gamma);
		SetXgamma (UtXgamma, UtX, rank_old);		
		logPost_old=CalcPosterior(UtXgamma, Utz, K_eval, UtXb_old, Utu_old, alpha_old, beta, cHyp_old);
	
		beta_old.clear();
		for (size_t i=0; i<beta->size; ++i) {
			beta_old.push_back(gsl_vector_get(beta, i));
		}	
		gsl_matrix_free (UtXgamma);
		gsl_vector_free (beta);
	}	
	
	//calculate centered z_hat, and pve
	if (a_mode==13) {
		time_start=clock();
		if (cHyp_old.n_gamma==0 || cHyp_old.rho==0) {
			CalcCC_PVEnZ (U, Utu_old, z_hat, cHyp_old);
		}
		else {
			CalcCC_PVEnZ (U, UtXb_old, Utu_old, z_hat, cHyp_old);
		}
		time_UtZ+=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);
	}
	
	//start MCMC
	int accept;
	size_t total_step=w_step+s_step;
	size_t w=0, w_col, pos;
	size_t repeat=0;
	
	for (size_t t=0; t<total_step; ++t) {
		if (t%d_pace==0 || t==total_step-1) {ProgressBar ("Running MCMC ", t, total_step-1, (double)n_accept/(double)(t*n_mh+1));}
//		if (t>10) {break;}		

		if (a_mode==13) {			
			SampleZ (y, z_hat, z);		
			mean_z=CenterVector (z);	
			
			time_start=clock();
			gsl_blas_dgemv (CblasTrans, 1.0, U, z, 0.0, Utz);
			time_UtZ+=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);
		
			//First proposal
			if (cHyp_old.n_gamma==0 || cHyp_old.rho==0) {				
				logPost_old=CalcPosterior(Utz, K_eval, Utu_old, alpha_old, cHyp_old);
				beta_old.clear();
				for (size_t i=0; i<cHyp_old.n_gamma; ++i) {
				  beta_old.push_back(0);
				}	
			}
			else {
				gsl_matrix *UtXgamma=gsl_matrix_alloc (ni_test, cHyp_old.n_gamma);
				gsl_vector *beta=gsl_vector_alloc (cHyp_old.n_gamma);
				SetXgamma (UtXgamma, UtX, rank_old);
				logPost_old=CalcPosterior(UtXgamma, Utz, K_eval, UtXb_old, Utu_old, alpha_old, beta, cHyp_old);
				
				beta_old.clear();
				for (size_t i=0; i<beta->size; ++i) {
					beta_old.push_back(gsl_vector_get(beta, i));
				}
				gsl_matrix_free (UtXgamma);
				gsl_vector_free (beta);
			}
		}
		
		//MH steps
		for (size_t i=0; i<n_mh; ++i) {
			if (gsl_rng_uniform(gsl_r)<0.33) {repeat = 1+gsl_rng_uniform_int(gsl_r, 20);}
			else {repeat=1;}
			
			logMHratio=0.0;
			logMHratio+=ProposeHnRho(cHyp_old, cHyp_new, repeat);		
			logMHratio+=ProposeGamma (rank_old, rank_new, p_gamma, cHyp_old, cHyp_new, repeat);	
			logMHratio+=ProposePi(cHyp_old, cHyp_new, repeat);
			
//			if (fix_sigma>=0) {
//				cHyp_new.h=fix_sigma/(1-cHyp_new.rho);
//			}
			
			if (cHyp_new.n_gamma==0 || cHyp_new.rho==0) {
				logPost_new=CalcPosterior(Utz, K_eval, Utu_new, alpha_new, cHyp_new);
				beta_new.clear();
				for (size_t i=0; i<cHyp_new.n_gamma; ++i) {
				  beta_new.push_back(0);
				}	
			}
			else {
				gsl_matrix *UtXgamma=gsl_matrix_alloc (ni_test, cHyp_new.n_gamma);
				gsl_vector *beta=gsl_vector_alloc (cHyp_new.n_gamma);
				SetXgamma (UtXgamma, UtX, rank_new);
				logPost_new=CalcPosterior(UtXgamma, Utz, K_eval, UtXb_new, Utu_new, alpha_new, beta, cHyp_new);
				beta_new.clear();
				for (size_t i=0; i<beta->size; ++i) {
					beta_new.push_back(gsl_vector_get(beta, i));
				}
				gsl_matrix_free (UtXgamma);
				gsl_vector_free (beta);
			}	
			
			logMHratio+=logPost_new-logPost_old;		
		
			if (logMHratio>0 || log(gsl_rng_uniform(gsl_r))<logMHratio) {accept=1; n_accept++;}
			else {accept=0;}
			
			if (accept==1) {			
				logPost_old=logPost_new;
				rank_old.clear(); beta_old.clear();
				if (rank_new.size()!=0) {
					for (size_t i=0; i<rank_new.size(); ++i) {
						rank_old.push_back(rank_new[i]);
						beta_old.push_back(beta_new[i]);
					}
				}
				cHyp_old=cHyp_new;
				gsl_vector_memcpy (alpha_old, alpha_new);
				gsl_vector_memcpy (UtXb_old, UtXb_new);
				gsl_vector_memcpy (Utu_old, Utu_new);
			}
			else {cHyp_new=cHyp_old;}
		}				
		
		//calculate z_hat, and pve
		if (a_mode==13) {
			time_start=clock();
			if (cHyp_old.n_gamma==0 || cHyp_old.rho==0) {
				CalcCC_PVEnZ (U, Utu_old, z_hat, cHyp_old);
			}
			else {
				CalcCC_PVEnZ (U, UtXb_old, Utu_old, z_hat, cHyp_old);
			}
			
			//sample mu and update z hat
			gsl_vector_sub (z, z_hat);
			mean_z+=CenterVector(z);
			mean_z+=gsl_ran_gaussian(gsl_r, sqrt(1.0/(double) ni_test) );			
			
			gsl_vector_add_constant (z_hat, mean_z);
			
			time_UtZ+=(clock()-time_start)/(double(CLOCKS_PER_SEC)*60.0);
		}
		
		//Save data
		if (t<w_step) {continue;}
		else {		
			if (t%r_pace==0) {
				w_col=w%w_pace;
				if (w_col==0) {
					if (w==0) {WriteResult (0, Result_hyp, Result_gamma, w_col);}					
					else {
						WriteResult (1, Result_hyp, Result_gamma, w_col);
						gsl_matrix_set_zero (Result_hyp);
						gsl_matrix_set_zero (Result_gamma);
					}
				}
				
				gsl_matrix_set (Result_hyp, w_col, 0, cHyp_old.h);
				gsl_matrix_set (Result_hyp, w_col, 1, cHyp_old.pve);
				gsl_matrix_set (Result_hyp, w_col, 2, cHyp_old.rho);
				gsl_matrix_set (Result_hyp, w_col, 3, cHyp_old.pge);
				gsl_matrix_set (Result_hyp, w_col, 4, cHyp_old.logp);
				gsl_matrix_set (Result_hyp, w_col, 5, cHyp_old.n_gamma);
				
				for (size_t i=0; i<cHyp_old.n_gamma; ++i) {
					pos=mapRank2pos[rank_old[i]]+1;

					gsl_matrix_set (Result_gamma, w_col, i, pos);
					
					beta_g[pos-1].first+=beta_old[i];
					beta_g[pos-1].second+=1.0;	
				}
				
				gsl_vector_add (alpha_prime, alpha_old);
				gsl_vector_add (Utu, Utu_old);
				
				if (a_mode==13) {
					pheno_mean+=mean_z;
				}
				
				w++;
				
			}
			
		}
	}
	cout<<endl;
	
	w_col=w%w_pace;
	WriteResult (1, Result_hyp, Result_gamma, w_col);	
	
	gsl_matrix_free(Result_hyp);
	gsl_matrix_free(Result_gamma);	
	
	gsl_vector_free(z_hat);
	gsl_vector_free(z);
	gsl_vector_free(Utz);	
	gsl_vector_free(UtXb_new);	
	gsl_vector_free(UtXb_old);
	gsl_vector_free(alpha_new);	
	gsl_vector_free(alpha_old);
	gsl_vector_free(Utu_new);	
	gsl_vector_free(Utu_old);	
	
	gsl_vector_scale (alpha_prime, 1.0/(double)w);	
	gsl_vector_scale (Utu, 1.0/(double)w);	
	if (a_mode==13) {
		pheno_mean/=(double)w;
	}
	
	gsl_vector *alpha=gsl_vector_alloc (ns_test);
	gsl_blas_dgemv (CblasTrans, 1.0/(double)ns_test, UtX, alpha_prime, 0.0, alpha);	
	WriteParam (beta_g, alpha, w);
	gsl_vector_free(alpha);
	
	gsl_blas_dgemv (CblasNoTrans, 1.0, U, Utu, 0.0, alpha_prime);
	WriteBV(alpha_prime);	
	
	gsl_vector_free(alpha_prime);				
	gsl_vector_free(Utu);	
		
	delete [] p_gamma;
	beta_g.clear();
	
	return;
}



void BSLMM::RidgeR(const gsl_matrix *U, const gsl_matrix *UtX, const gsl_vector *Uty, const gsl_vector *eval, const double lambda)
{
	gsl_vector *beta=gsl_vector_alloc (UtX->size2);
	gsl_vector *H_eval=gsl_vector_alloc (Uty->size);
	gsl_vector *bv=gsl_vector_alloc (Uty->size);

	gsl_vector_memcpy (H_eval, eval);
	gsl_vector_scale (H_eval, lambda);
	gsl_vector_add_constant (H_eval, 1.0);
	
	gsl_vector_memcpy (bv, Uty);
	gsl_vector_div (bv, H_eval);	

	gsl_blas_dgemv (CblasTrans, lambda/(double)UtX->size2, UtX, bv, 0.0, beta);
	gsl_vector_add_constant (H_eval, -1.0);
	gsl_vector_mul (H_eval, bv);
	gsl_blas_dgemv (CblasNoTrans, 1.0, U, H_eval, 0.0, bv);

	WriteParam (beta);
	WriteBV(bv);
	
	gsl_vector_free (H_eval);
	gsl_vector_free (beta);
	gsl_vector_free (bv);
	
	return;
}
 
