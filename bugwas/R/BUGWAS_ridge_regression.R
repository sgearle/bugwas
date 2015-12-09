# File KmerAnalysisMain.R
# Authors: Earle, S. G., Wu, C.-H. and Wilson, D. J.
#
# Copyright (C) 2015 University of Oxford
#
# This file is part of the bugwas R package.
# See the NOTICE file distributed with this work for additional
# information regarding copyright ownership and licensing.
#
# This is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
#  This is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this software package; if not, write to the
# Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
# Boston, MA  02110-1301  USA

ridge_regression = function(y,x,w=matrix(1,length(y),1),K_type="uncentred",svdX=NULL,lambda_init=200,prefer.REML=TRUE,skip.var=FALSE,maximize=TRUE) {
	if(lambda_init<=0) stop("ridge_regression: lambda_init must be positive")
	n = length(y)
# Process x as uncentred, centred or standardized, and label the result X
# (not yet implemented: take x at face value)
	if(K_type!="uncentred") stop("Currently the data x must be pre- centred or standardized if desired")
	X = x
#	Note that after centring/standardizing X, K = X t(X) (i.e. K_ij = sum_l X_il X_jl)
#	and that by singular value decomposition (SVD), X = U D V`, with U and V orthogonal,
#	i.e. V^(-1) = V`. So one can write
#	K = X t(X) = U D V` V` D U` = U D^2 U`.
#	If you define R = U D, then K = X t(X) = R t(R)
#	Calculating K via SVD will be faster in general because R is only an n x n matrix
# Do the singular value decomposition if necessary
	if(is.null(svdX)) {
		svdX = svd(X)
	}
	rk = length(svdX$d)
	if(rk<n) {
		svdX$d = c(svdX$d,rep(0,n-rk))
		svdX$u = cbind(svdX$u,matrix(0,n,n-rk))
		svdX$v = cbind(svdX$v,matrix(0,nrow(svdX$v),n-rk))
	}
# Calculate the relatedness matrix
	svdR = svdX$u %*% diag(svdX$d)
#	svdR = matrix(0,n,n)
#	svdR[1:n,1:length(svdX$d)] = svdX$u %*% diag(svdX$d)
	K = (svdR %*% t(svdR)) # Unlike GEMMA, not dividing by number of loci

# Precompute other quantities
	Wx = as.matrix(w)
	c = ncol(w)

# Estimate lambda using maximum likelihood
	lambda_loglik = function(log_lambda) {
		lambda = exp(log_lambda)
		H = lambda * K + diag(n)
		Hinv = solve(H)
		Px = Hinv - (Hinv %*% Wx %*% solve(t(Wx) %*% Hinv %*% Wx) %*% t(Wx) %*% Hinv)
		ret = 0.5*n*log(0.5*n/pi) - 0.5*n - 0.5*as.numeric(determinant(H,log=TRUE)$modulus) - 0.5*n*log(t(y) %*% Px %*% y)
		if(!is.finite(ret)) return(-.Machine$double.xmax)
		return(ret)
	}
#(lambda_MLE = optimize(lambda_loglik,c(1e-5,1e5),maximum=TRUE)$maximum)
#	nlm_MLE = suppressWarnings(nlm(function(x) -lambda_loglik(x),log(lambda_init)))
	if(maximize) {
		nlm_MLE = suppressWarnings(nlm(function(x) tryCatch(-lambda_loglik(x),error=function(e) -Inf),log(lambda_init)))
	} else {
		nlm_MLE = list("estimate"=log(lambda_init),"minimum"=-lambda_loglik(log(lambda_init)))
	}
	lambda_MLE = exp(nlm_MLE$estimate)
	ML = -nlm_MLE$minimum
# Estimate lambda using REML
	lambda_logpartiallik = function(log_lambda) {
		lambda = exp(log_lambda)
		H = lambda * K + diag(n)
		Hinv = solve(H)
		Px = Hinv - (Hinv %*% Wx %*% solve(t(Wx) %*% Hinv %*% Wx) %*% t(Wx) %*% Hinv)
# Note that whereas Zhou and Stephens use (n-c-1) we use (n-c) because we do not include an extra fixed effect for an individual genotype
		ret = 0.5*(n-c)*log(0.5*(n-c)/pi) - 0.5*(n-c) + 0.5*as.numeric(determinant(t(Wx)%*%Wx,log=TRUE)$modulus) - 0.5*as.numeric(determinant(H,log=TRUE)$modulus) - 0.5*as.numeric(determinant(t(Wx)%*%Hinv%*%Wx,log=TRUE)$modulus) - 0.5*(n-c)*log(t(y) %*% Px %*% y)
		if(!is.finite(ret)) return(-.Machine$double.xmax)
		return(ret)
	}
#(lambda_REML = optimize(lambda_logpartiallik,c(1e-5,1e5),maximum=TRUE)$maximum)
#	nlm_REML = suppressWarnings(nlm(function(x) -lambda_logpartiallik(x),log(lambda_init)))
	if(maximize) {
		nlm_REML = suppressWarnings(nlm(function(x) tryCatch(-lambda_logpartiallik(x),error=function(e) -Inf),log(lambda_init)))
	} else {
		nlm_REML = list("estimate"=log(lambda_init),"minimum"=-lambda_logpartiallik(log(lambda_init)))
	}
	lambda_REML = exp(nlm_REML$estimate)
	REML = -nlm_REML$minimum
	
# Estimate a and tau given lambda using maximum likelihood
	lambda = lambda_MLE
	H = lambda * K + diag(n)
	Hinv = solve(H)
	a_MLE = solve(t(Wx) %*% Hinv %*% Wx) %*% t(Wx) %*% Hinv %*% y
	Px = Hinv - (Hinv %*% Wx %*% solve(t(Wx) %*% Hinv %*% Wx) %*% t(Wx) %*% Hinv)
	tau_MLE = n / (t(y) %*% Px %*% y)
# Estimate tau given lambda using REML
	lambda = lambda_REML
	H = lambda * K + diag(n)
	Hinv = solve(H)
	a_MLE = solve(t(Wx) %*% Hinv %*% Wx) %*% t(Wx) %*% Hinv %*% y
	Px = Hinv - (Hinv %*% Wx %*% solve(t(Wx) %*% Hinv %*% Wx) %*% t(Wx) %*% Hinv)
# Note that whereas Zhou and Stephens use (n-c-1) we use (n-c) because we do not include an extra fixed effect for an individual genotype
	tau_REML = (n - c) / (t(y) %*% Px %*% y)
	
# Calculate the posterior mean and variance in effect size for each SNP
	lambda = lambda_REML
	tau = tau_REML
	if(!prefer.REML) {
		lambda = lambda_MLE
		tau = tau_MLE
	}
# Note that with prior mean m=0, the posterior mean is the ridge regression estimate (Eq 11.57 of O`Hagan & Forster with tau=sigma^(-2))
#	mstar	= solve(t(X) X + 1/lambda I_L) t(X) y
# where L is the number of loci. This would involve inverting a very large L x L matrix. It is more efficient to use SVD:
#	mstar	= solve(X` X + 1/lambda I_L) X` y
# next introduce identity matrices
#			= (V solve(V)) solve(X` X + 1/lambda I_L) (solve(V`) V`) X` y
# then use the rule solve(A) solve(B) solve(C) = solve(C B A)
#			= V solve(V` (X` X + 1/lambda I_L) V) V` X` y
# next expand the brackets
#			= V solve(V` X` X V + 1/lambda V` I_L V) V` X` y
# and note that V` I_L V = I_n
#			= V solve(V` X` X V + 1/lambda I_n) V` X` y
# using SVD, write X = U D V` and X` = V D U`
#			= V solve(V` V D U` U D V` V + 1/lambda I_n) V` V D U` y
# and use the property that V` V = I_n
#			= V solve(D U` U D + 1/lambda I_n) D U` y
# then introduce R = U D and R` = D U`
#			= V solve(R` R + 1/lambda I_n) R` y
# finally define C = solve(R` R + 1/lambda I_n) to give
#			= V C R` y
	C = solve(t(svdR) %*% svdR + (1/lambda * diag(n)))
	post.mean = (svdX$v %*% (C %*% t(svdR) %*% y))
# The posterior variance, Wstar = solve(solve(W) + tau t(X) X) (Eq 11.61 of O`Hagan & Forster, p327 with tau=sigma^(-2))
# Assuming a prior variance W = lambda/tau I_L,
#	Wstar	= lambda/tau solve(I_L + lambda t(X) X)
# Applying the SVD,
#			= lambda/tau solve(I_L + lambda V D U` U D V`)
#			= lambda/tau solve(I_L + lambda V D D V`)
# substituting E = lambda D^2 (a diagonal matrix)
#			= lambda/tau solve(I_L + lambda V E V`)
# Following section 1.5 of Henderson and Searle 1981 SIAM Review 23:53-60, Eq. 17, with A = I_n, B = E and U = V. In their notation,
#			solve(A + U B U`) = solve(A) - solve(A) U solve(solve(B) + U` solve(A) U) U` solve(A)
# gives in our notation
#			solve(I_L + V E V`) = I_L - I_L V solve(solve(E) + V` I_L V) V` I_L = I_L - V solve(solve(E) + I_n) V`
# so
#	Wstar	= lambda/tau (I_L - V Cstar V`)
# where
#	Cstar	= lambda D^2/(lambda D^2 + 1) I_n
	Cstar = diag(lambda * svdX$d^2 / (lambda * svdX$d^2 + 1))
	if(skip.var) {
		post.var = NA 
	} else {
#		post.var = lambda/tau * apply(svdX$v,1,function(V) 1 - V %*% Cstar %*% V)
#		post.var = lambda/tau * sapply(1:nrow(svdX$v),function(i) as.vector(1 - svdX$v[i,] %*% Cstar %*% svdX$v[i,]))
#		post.var = lambda/tau * sapply(1:nrow(svdX$v),function(i) 1 - sum(svdX$v[i,]^2*diag(Cstar)))
		post.var = lambda/tau * (1 - colSums(t(svdX$v^2)*diag(Cstar)))
	}
	
	ret = list("lambda_MLE"=lambda_MLE,"lambda_REML"=lambda_REML,"a_MLE"=a_MLE,"tau_MLE"=tau_MLE[1,1],"tau_REML"=tau_REML[1,1],"Ebeta"=post.mean,"Vbeta"=post.var,"ML"=ML,"REML"=REML,"func"=function(x){x})
	return(ret)
}