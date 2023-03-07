#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h>
#include <vector>
#include <string>
#include <random>
#include <stdlib.h>
// [[Rcpp::depends(BH)]]

#include <boost/random.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/math/distributions.hpp>
#include <boost/random/variate_generator.hpp>
#include <algorithm>
// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace std;
using namespace arma;
using namespace Rcpp;
boost::random::mt19937 rng;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]

mat submatrix(mat y, uvec v){
  uvec idx = find(v == 1);
  return y.cols(idx);
}


arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}


double rInvGamma(const double shape, const double scale){
  boost::gamma_distribution<> gd( shape );
  boost::variate_generator<boost::mt19937&, boost::gamma_distribution<> > var_gamma( rng, gd );
  double out = 0;
  out = scale/(var_gamma());
  //we don't use randg because it doesn't take in double shape and scale.
  return out;
}


double rExp(const double lbda ){
  double randFromUnif = randu();
  //parameters and the random value on (0,1) you drew
  boost::math::exponential_distribution<> dist(lbda);
  double randFromDist = quantile(dist, randFromUnif);
  return randFromDist;
}


double dExp(const double x, const double lbda ){
  boost::math::exponential_distribution<> my_exp(lbda);
  double out = boost::math::pdf(my_exp,x);
  return out;
}


// proposal of different gamma
uvec proposal(uvec gamma, double rho){
  int k = sum(gamma);
  int len = gamma.size();
  if(len == 1){
    return 1-gamma;
  }
  uvec gamma_new = repelem(gamma,1,1);
  if(k == 0){
    int idx = as_scalar(randperm((len-1),1));
    
    gamma_new(idx) = 1;
  }
  else if(k == len){
    int idx = as_scalar(randperm((len-1),1));
    
    gamma_new(idx) = 0;
  }
  else{
    double tmpU = as_scalar(randu());
    if(tmpU < rho){
      int idx = randi(distr_param(0,len-1));
      int v = 0.5 + 0.5*pow(-1, gamma(idx));
      gamma_new(idx) = v;
    }
    else{
      uvec gamma1 = find(gamma == 1);
      uvec gamma2 = find(gamma == 0);
      uvec gamma1_perm = gamma1(randperm(gamma1.size()));
      uvec gamma2_perm = gamma2(randperm(gamma2.size()));
      gamma_new(gamma1_perm(0)) = 0;
      gamma_new(gamma2_perm(0)) = 1;
    }
  }
  return gamma_new;
}


// [[Rcpp::export]]
double clinic_loglikelihood_integrate(vec Y, mat X, uvec gamma, mat Hbar, double tau, double delta1, double delta2){
  
  int n = Y.size();
  int k = sum(gamma);
  double Q = as_scalar(Y.t() * Hbar * Y);
  
  if(k == 0){
    
    double out = -(0.5 * n + delta1) * log(0.5 * Q + delta2);
    return out;
  }
  else{
    double term1 = 0.5 * k * log(tau);
    mat X_tuple = submatrix(X, gamma);
    mat A;
    mat H_tuple = X_tuple.t() * Hbar * X_tuple + tau*A.eye(k,k); 
    double term2 = - 0.5 * log_det(H_tuple).real();
    mat P = solve(H_tuple, X_tuple.t() * Hbar * Y);
    double Q2 = as_scalar(Y.t() * Hbar * X_tuple * P);
    Q = Q -  Q2 + 0.0001; 
    double term3 = -(0.5*n + delta1)*log(0.5*Q + delta2);
    return term1 + term2 + term3;
  }
}

// [[Rcpp::export]]
double conditional_loglikelihood_integrate(vec B, double beta1, double beta2, vec g, mat M, uvec Zj, double sigma2, double sigma_g2, double tau){
  int k = sum(Zj);
  if(k == 0)
  {
    return(0);
  }
  else{
    double term1 = 0.5 * k * log(tau) - 0.5*log(sigma_g2);
    mat M_tuple = submatrix(M, Zj);
    mat M_cov = M_tuple.t() * M_tuple;
    mat A;
    mat pres = 1/sigma2 * pow((beta2 - beta1),2) * M_cov + 1/sigma_g2 * (M_cov + tau*A.eye(k,k));
    double term2 = - 0.5 * log_det_sympd(pres);
    mat X = 1/sigma2 *(beta2 - beta1) * M_tuple.t() * B + 1/sigma_g2 * M_tuple.t() * g;
    mat P = solve(pres, X);
    double term3 = as_scalar(0.5 * X.t() * P);
    return term1 + term2 +term3;
  }
}

// [[Rcpp::export]]
double log_MRF(uvec gamma1, uvec gamma2, double a, double b, mat G){
  uvec gamma_u = gamma1 + gamma2;
  double out = a * sum(gamma_u);
  gamma_u.elem(find(gamma_u > 0)).ones();
  out = out + b * as_scalar(gamma_u.t() * G *gamma_u);
  return out;
}


// Gibbs algorithm to update regression coefficients of genes
// [[Rcpp::export]]
vec Gibbs_update_beta(vec Y, mat X, uvec gamma, double tau, double sigma2){

  mat Z = submatrix(X, gamma);
  int k = sum(gamma);
  vec B = tau * ones<vec>(k);
  mat A = diagmat(B);
  mat precision = Z.t() * Z + A;
  mat covariance = inv_sympd(precision);
  vec mu = vectorise(covariance * Z.t() * Y);
  mat beta = mvrnormArma(1, mu, covariance*sigma2);
  return vectorise(beta);
}

// Gibbs algorithm to update regression coefficients of clinic covariates
// [[Rcpp::export]]
vec Gibbs_update_beta_c(vec Y, mat C, double tau_c, double sigma2){
  
  int k = C.n_cols;
  vec B = tau_c * ones<vec>(k);
  mat A = diagmat(B);
  mat precision = C.t() * C + A;
  mat covariance = inv_sympd(precision);
  vec mu = vectorise(covariance * C.t() * Y);
  mat beta = mvrnormArma(1, mu, covariance*sigma2);
  return vectorise(beta);
}


// Gibbs algorithm to update Omega_j
// [[Rcpp::export]]
vec Gibbs_update_omega(vec B, double beta1, double beta2, vec g, mat M, uvec Zj, double sigma2, double sigma_g2, double tau){

  int k = sum(Zj);
  int p = Zj.size();
  if(k == 0 )
  {
    return zeros<vec>(p);
  }
  mat M_tuple = submatrix(M, Zj);
  mat M_cov = M_tuple.t() * M_tuple;
  mat pres = 1/sigma2 * pow((beta2 - beta1),2) * M_cov + 1/sigma_g2 * (M_cov + tau*eye<mat>(k,k));
  mat X = 1/sigma2 *(beta2 - beta1) * M_tuple.t() * B + 1/sigma_g2 * M_tuple.t() * g;
  vec mu = vectorise(solve(pres, X));
  mat var = inv_sympd(pres);
  vec omega = vectorise(mvrnormArma(1, mu, var));
    return omega;
}

// Gibbs algorithm to update sigma2
// [[Rcpp::export]]
double Gibbs_update_sigma2(vec Y, mat X, vec beta, double delta1, double delta2)
{
  int n = Y.size();
  vec err = vectorise(Y - X * beta);
  err = pow(err, 2);
  double result = rInvGamma((delta1 + n/2.0), (delta2 + sum(err)/2.0));
  return result;
}


// MH algorithm to update gamma
// [[Rcpp::export]]
uvec MH_update_gamma(vec Y, mat X, uvec gamma, mat Hbar, double tau, double delta1, double delta2, double a, double b, mat graph, double rho){

  int p = graph.n_cols;
  uvec gamma_new = join_cols(proposal(gamma.subvec(0, (p-1)), rho), proposal(gamma.subvec(p, (2*p-1)), rho));

  double Term1 = clinic_loglikelihood_integrate(Y, X, gamma_new, Hbar, tau, delta1, delta2) - clinic_loglikelihood_integrate(Y, X, gamma, Hbar, tau, delta1, delta2);
  double Term2 = log_MRF(gamma_new.subvec(0, (p-1)), gamma_new.subvec(p, (2*p-1)), a, b, graph) - log_MRF(gamma.subvec(0, (p-1)), gamma.subvec(p, (2*p-1)), a, b, graph);
  double tmpU = as_scalar(randu());
  double thres = Term1 + Term2;
  vec temp = {thres, -999.0};
  if(temp.max() > log(tmpU))
  {
    return gamma_new;
  }
  else{return gamma;}
}



// Gibbs algorithm to update gamma
// [[Rcpp::export]]
uvec Gibbs_update_gamma(vec R, mat X, uvec gamma, vec beta, double sigma2, double tau, double a, double b, mat graph)
{
  int p = graph.n_cols;
  uvec gamma_new = zeros<uvec>(2*p);
  
  for(int j=0;j<(2*p);j++){
    vec res = R + X.col(j) * beta(j);
    vec temp = pow(X.col(j), 2);
    double factor1 = sum(temp);
    double factor2 = pow(sum(res % X.col(j)), 2);
    double log_ratio = 0.5*(log(tau) - log(tau + factor1)) + factor2 / (sigma2*(tau + factor1));
    uvec gamma1 = repelem(gamma,1,1);
    gamma1(j) = 1;
    double Term2 = log_MRF(gamma1.rows(0,(p-1)), gamma1.rows(p,(2*p-1)), a, b, graph) - log_MRF(gamma.rows(0,(p-1)), gamma.rows(p,(2*p-1)), a, b, graph);
    double prob = min(exp(log_ratio)*exp(Term2), 999.0);
    prob = prob / (1.0 + prob);
    double tmpU = as_scalar(randu());
    if(tmpU < prob){
      gamma_new(j) = 1;
    }
    else{
      gamma_new(j) = 0;
    }
  }
  return gamma_new;
}


// MH algorithm to update Zj
// [[Rcpp::export]]
uvec MH_update_Z(vec B, double beta1, double beta2, vec g, mat M, uvec Zj, double sigma2, double sigma_g2, double tau, double e, double f, double rho){

  int p = Zj.size();
  uvec Zj_new = proposal(Zj, rho);
  double Term1 = conditional_loglikelihood_integrate(B, beta1, beta2, g, M, Zj_new, sigma2, sigma_g2, tau) - conditional_loglikelihood_integrate(B, beta1, beta2, g, M, Zj, sigma2, sigma_g2, tau);
  double Term2 = 0;
  for(int i=0;i< p;i++){
    Term2 += log(tgamma(e + Zj_new(i))) + log(f +1.0 - Zj_new(i)) - log(tgamma(e + Zj(i))) -log(tgamma(f + 1.0 -Zj(i)));
  }
  double tmpU = as_scalar(randu());
  double thres = Term1 + Term2;
  vec temp = {thres, -999.0};
  if(temp.max() > log(tmpU))
  {
    return Zj_new;
  }
  else{return Zj;}
}


// missing gene expression sampling algorithm 
// [[Rcpp::export]]
vec mRNA_update(vec A, mat M, double beta, vec Omega, double sigma2, double sigma_g2)
{
  
  int n_miss = A.size();
  double variance = 1/(pow(beta, 2)/sigma2 + sigma_g2);
  vec mu = variance * vectorise( (A*beta / sigma2 + M * Omega /sigma_g2));
  mat temp;
  vec output = vectorise(mvrnormArma(1, mu, variance * temp.eye(n_miss, n_miss)));
  return output;
}

// missing DNA methylation sampling algorithm
// [[Rcpp::export]]
mat methyl_update(vec A, vec E, mat M, vec Omega, double beta1, double beta2, double sigma2, double sigma_g2)
{
  int q = Omega.size();
  int n_miss = A.size();
  mat B = E - M * Omega;
  mat output = zeros<mat>(n_miss, q);
  for(int i=0;i<q;i++){
    
    double precision = 1.0 + (pow((beta2 - beta1),2) / sigma2 + 1/ sigma_g2)*pow(Omega(i),2);
    vec mu = 1.0/precision * (A * (beta2-beta1)/sigma2 + (B + M.col(i)*Omega(i))/sigma_g2)*Omega(i);
    mat temp;
    output.col(i) = vectorise(mvrnormArma(1, mu,  1.0/precision * temp.eye(n_miss, n_miss)));
  }
  return output;
}


// [[Rcpp::export]]
mat methyl_block_product(mat M, Rcpp::List& list_map, vec Omega){
  int p = list_map.size();
  int n = M.n_rows;
  mat out = zeros<mat>(n, p);
  for(int j=0; j<p; j++){
    uvec idx = Rcpp::as<uvec>(list_map(j)) -1 ;
    out.col(j) = M.cols(idx) * Omega.rows(idx);
  }
  return out;
}


// [[Rcpp::export]]
double MH_update_Y(int index, vec Y_censor, uvec censorindex, vec Y, mat X, uvec gamma, mat Hbar, double tau, double delta1, double delta2){
  
  int pos = censorindex(index-1);
  double proposal = rExp(1/(Y(pos-1) - Y_censor(index-1))) + Y_censor(index-1);
  vec Y_new = Y;
  Y_new(pos-1) = proposal;
  
  double Term1 = clinic_loglikelihood_integrate(Y_new, X, gamma, Hbar, tau, delta1, delta2) - clinic_loglikelihood_integrate(Y, X, gamma, Hbar, tau, delta1, delta2);
  double Term2 = log(dExp(Y(pos-1)- Y_censor(index-1), 1/(Y_new(pos-1) - Y_censor(index-1)))) - log(dExp(Y_new(pos-1)- Y_censor(index-1), 1/(Y(pos-1) - Y_censor(index-1))));
  
  double tmpU = as_scalar(randu());
  double thres = Term1 + Term2;
  vec temp = {thres, -999.0};
  if(temp.max() > log(tmpU))
  {
    return proposal;
  }
  else{return Y(pos-1);}
  
}



