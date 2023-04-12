# call utils.cpp function
# final version applicable for missing data

library(Rcpp)
library(RcppArmadillo)
library(RcppGSL)
library(scales)
library(MASS)


Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp("R/utils.cpp")

### posterior of clinical model
posterior_clinic_model <- function(Y, C, X, beta_c, beta, gamma, sigma2, tau_c, tau, a, b, graph, delta1, delta2)
{
  n = length(Y)
  p = as.integer(length(gamma)/2)
  k = sum(gamma)
  Term1 = -0.5 * n * log(sigma2) - 0.5 * sum((Y - C %*% beta_c - X %*% beta)^2)/sigma2
  Term2 = -0.5 * k * log(2*pi) - 0.5 * k * log(sigma2) + 0.5 * k * log(tau) - 0.5 * tau * sum(beta^2)/sigma2
  Term3 = log_MRF(gamma[1:p], gamma[-(1:p)], a, b, graph)
  Term4 = -0.5 * tau_c * sum(beta_c^2)/sigma2
  Term5 = -delta1*log(sigma2) - delta2/sigma2
  return(Term1+Term2+Term3+Term4+Term5)
}

### likelihood of clinical model
likelihood_clinic_model <- function(Y, C, X, beta_c, beta, sigma2)
{
  n = length(Y)
  Term1 = -0.5 * n * log(sigma2) - 0.5 * sum((Y - C %*% beta_c - X %*% beta)^2)/sigma2
  return(Term1)
}

### log density for each component of  Zj
log.density.z <- function(z, e, f){
  return(log(gamma(e+z)) + log(gamma(f + 1 -z)))
}


### posterior of mechanistic model
posterior_mechanistic_model <- function(E, M, Omega, Zj, sigma_g2, tau_g, e, f, delta1, delta2)
{
  Omega = as.matrix(Omega)
  k = sum(Zj)
  n = length(E)
  Term1 =  - 0.5 * n * log(sigma_g2) - 0.5 * sum((E - M %*% Omega )^2)/sigma_g2
  Term2 = -0.5 * k * log(2*pi) - 0.5 * k * log(sigma_g2) + 0.5 * k * log(tau_g) - 0.5 * tau_g * sum(Omega^2)/sigma_g2
  Term3 = -delta1*log(sigma_g2) - delta2/sigma_g2
  Term4 = sum(log.density.z(Zj, e, f))
  return(Term1 + Term2 + Term3 + Term4)
}


### likelihood of mechanistic model
likelihood_mechanistic_model <- function(E, M, Omega, sigma_g2)
{
  n = length(E)
  Omega = as.matrix(Omega)
  Term1 =  - 0.5 * n * log(sigma_g2) - 0.5 * sum((E - M %*% Omega )^2)/sigma_g2
  return(Term1)
}




NIBAMM <- function(data, init=NULL, prior, max_iters, seed=15267){
  set.seed(seed)
  gammaTimes = 10
  rho = 0.5
  # extract data
  Y = data$Y
  C = as.matrix(data$C)
  E = as.matrix(data$E)
  M = as.matrix(data$M)
  map = data$map
  
  p1 = dim(C)[2]
  p = dim(E)[2]
  n = length(Y)
  p2 = dim(M)[2]
  
  list_map = list()
  for(j in 1:p){
    list_map = append(list_map, list(which(map == j)))
  }
  
  # hyper parameter in the prior
  tau_c = prior$tau_c
  tau = prior$tau
  tau_g = prior$tau_g
  delta1 = prior$delta1
  delta2 = prior$delta2
  a = prior$a
  b = prior$b
  e = prior$e
  f = prior$f
  graph = prior$graph
  
  ## data structure to store parameter
  beta_c = matrix(0, nrow = max_iters, ncol = p1)
  beta = matrix(0, nrow = max_iters, ncol = 2*p)
  gamma = matrix(NA, nrow = max_iters*gammaTimes, ncol = 2*p)
  Omega = matrix(0, nrow = max_iters, ncol = p2)
  Z = matrix(NA, nrow = max_iters, ncol = p2)
  sigma2 = rep(NA, max_iters)
  sigma_g2 = matrix(NA, nrow = max_iters, ncol = p)
  posterior = rep(0, max_iters)
  likelihood = rep(0, max_iters)
  
  # determine missing samples
  miss_exp_sample = which(rowSums(matrix(as.numeric(is.na(E)), n, p)) > 0)
  n_exp_miss = length(miss_exp_sample)
  
  miss_met_sample = which(rowSums(matrix(as.numeric(is.na(M)), n, p2)) > 0)
  n_met_miss = length(miss_met_sample)
  
  if(is.null(init)){
    # initial of the parameter
    beta_c[1, ] = rnorm(p1)
    beta[1, ] = rnorm(2*p)
    gamma[1:gammaTimes, ] = matrix(rep(0, gammaTimes*2*p), gammaTimes, 2*p)
    Omega[1, ] = rnorm(p2)
    Z[1, ] = rep(0, p2)
    sigma2[1] = rgamma(1, shape = 1, rate = 1)
    sigma_g2[1, ] = rgamma(p, shape = 1, rate = 1)
    posterior[1] = -Inf
    likelihood[1] = -Inf
  }
  else{
    # initial of the parameter
    beta_c[1, ] = init$beta_c
    beta[1, ] = init$beta
    gamma[1:gammaTimes, ] = matrix(rep(init$gamma, gammaTimes), gammaTimes, 2*p)
    Omega[1, ] = init$Omega
    Z[1, ] = init$Z
    sigma2[1] = init$sigma2
    sigma_g2[1, ] = init$sigma_g2
    posterior[1] = -Inf
    likelihood[1] = -Inf
  }
  ## initialize missing data
  if(n_exp_miss != 0){
    
    ### initial of missing expression data
    #E[miss_exp_sample, ] = matrix(rep(colMeans(E, na.rm = T), each=n_exp_miss), nrow = n_exp_miss, ncol = p)
    E[miss_exp_sample, ] = matrix(rnorm(n_exp_miss*p), nrow = n_exp_miss, ncol = p)
  }
  
  if(n_met_miss != 0){
    ### initial of missing methylation data
    M[miss_met_sample, ] = matrix(rep(colMeans(M, na.rm = T), each=n_met_miss), nrow = n_met_miss, ncol = p2)
  }
  
  # some important intermediate variable
  Hhat =  t(C) %*% C + diag(tau_c, nrow = p1)
  Hhat_inverse = solve(Hhat, t(C))
  Hbar = diag(1, n) - C %*% Hhat_inverse 
  
  Ehat = methyl_block_product(M, list_map, Omega[1, ])
  
  X = matrix(0, n ,2*p)
  
  if(n_exp_miss !=0 & n_met_miss !=0){
    
    # MCMC process
    for(i in 2:max_iters)
    {
      
      ### update parameter in mechanistic model
      B_tutle = Y - C %*% beta_c[i-1, ] - E %*% beta[i-1, 1:p]
      beta1 = beta[i-1, 1:p]
      beta2 = beta[i-1, -(1:p)]
      for(j in 1:p)
      {
        B = B_tutle - Ehat[, -j] %*% (beta2[-j] - beta1[-j])
        temp = as.matrix(M[, list_map[[j]]])
        
        Z[i, list_map[[j]]] = MH_update_Z(B=B, beta1[j], beta2[j], g=E[, j], temp, Zj=Z[i-1, list_map[[j]]], sigma2[i-1], sigma_g2[i-1, j], tau_g[j],e, f, rho=rho)
        
        if(sum(Z[i, list_map[[j]]]) != 0){
          updated_value = Gibbs_update_omega(B=B, beta1[j], beta2[j], g=E[, j], temp, Zj=Z[i, list_map[[j]]], sigma2[i-1], sigma_g2[i-1, j], tau_g[j])
          idx = list_map[[j]][Z[i, list_map[[j]]] == 1]
          Omega[i, idx] = updated_value
        }
        else{
          Omega[i, list_map[[j]]] = 0
        }
        
        sigma_g2[i, j] = Gibbs_update_sigma2(E[, j], temp, Omega[i,list_map[[j]]], delta1, delta2)
        
        posterior[i] = posterior[i] + posterior_mechanistic_model(E[,j], M[ , list_map[[j]]], Omega[i,list_map[[j]] ], Zj=Z[i, list_map[[j]]], sigma_g2[i, j], tau_g[j], e, f, delta1, delta2)
        likelihood[i] = likelihood[i] + likelihood_mechanistic_model(E[,j], M[ , list_map[[j]]], Omega[i,list_map[[j]] ], sigma_g2[i, j])
      }
      
      ### update parameter in clinical model
      Ehat = methyl_block_product(M, list_map, Omega = Omega[i, ])
      X = cbind(E-Ehat, Ehat)
      
      for(k in 1:gammaTimes){
        gamma[(i-1)*gammaTimes + k , ] = MH_update_gamma(Y, X, gamma[(i-1)*gammaTimes + k-1, ], Hbar, tau, delta1, delta2, a, b, graph, rho=rho)
        
      }
      
      ### update beta
      err = Y - C %*% beta_c[i-1, ]
      updated_value = Gibbs_update_beta(err, X, gamma[i*gammaTimes, ], tau, sigma2[i-1])
      
      beta[i, gamma[i*gammaTimes, ] == 1] = updated_value
      beta[i, gamma[i*gammaTimes, ] != 1] = 0
      
      ### update beta_c
      err = Y - X %*% beta[i, ]
      beta_c[i, ] = Gibbs_update_beta_c(err, C, tau_c, sigma2[i-1] )
      
      ### update variance
      sigma2[i] = Gibbs_update_sigma2(Y, cbind(C, X), c(beta_c[i, ], beta[i, ]), delta1, delta2)
      
      posterior[i] = posterior[i] + posterior_clinic_model(Y, C, X, beta_c[i,], beta[i,], gamma[i*gammaTimes, ], sigma2[i], tau_c, tau, a, b, graph, delta1, delta2)
      likelihood[i] = likelihood[i] + likelihood_clinic_model(Y, C, X, beta_c[i,], beta[i,], sigma2[i])
      
      ### update missing data
      E_mis_bar = E[miss_exp_sample, ] %*% beta[i, 1:p] 
      E_mis_hat = Ehat[miss_exp_sample, ] %*% (beta[i, -(1:p)] - beta[i, 1:p])
      A_mis = Y[miss_exp_sample] - C[miss_exp_sample, ] %*% beta_c[i, ] - E_mis_hat - E_mis_bar
      for(j in 1:p)
      {
        beta_j = beta[i,1:p][j]
        A = A_mis + E[miss_exp_sample,j] * beta_j
        
        E[miss_exp_sample, j] = mRNA_update(A=A, M=as.matrix(M[miss_exp_sample, list_map[[j]]]), beta_j, Omega[i, list_map[[j]] ], sigma2[i], sigma_g2[i,j])
      }
      
      E_mis_bar = E[miss_met_sample, ] %*% beta[i, 1:p] 
      E_mis_hat = Ehat[miss_met_sample, ] %*% (beta[i, -(1:p)] - beta[i, 1:p])
      A_mis = Y[miss_met_sample] - C[miss_met_sample, ] %*% beta_c[i, ] - E_mis_hat - E_mis_bar
      for(j in 1:p)
      {
        beta_j1 = beta[i, 1:p][j]
        beta_j2 = beta[i, -(1:p)][j]
        A = A_mis + as.matrix(M[miss_met_sample, list_map[[j]] ]) %*% Omega[i, list_map[[j]]]
        updated_value = methyl_update(A, E[miss_met_sample, j], as.matrix(M[miss_met_sample, list_map[[j]]]), Omega[i,list_map[[j]]], beta_j1, beta_j2, sigma2[i], sigma_g2[i, j])
        M[miss_met_sample, list_map[[j]]] = updated_value
      }
    }
  }
  else if(n_exp_miss ==0 & n_met_miss !=0){
    # MCMC process
    for(i in 2:max_iters)
    {
      ### update parameter in mechanistic model
      B_tutle = Y - C %*% beta_c[i-1, ] - E %*% beta[i-1, 1:p]
      beta1 = beta[i-1, 1:p]
      beta2 = beta[i-1, -(1:p)]
      for(j in 1:p)
      {
        B = B_tutle - Ehat[, -j] %*% (beta2[-j] - beta1[-j])
        temp = as.matrix(M[, list_map[[j]]])
        Z[i, list_map[[j]]] = MH_update_Z(B=B, beta1[j], beta2[j], g=E[, j], temp, Zj=Z[i-1, list_map[[j]]], sigma2[i-1], sigma_g2[i-1, j], tau_g[j],e, f, rho=rho)
        
        if(sum(Z[i, list_map[[j]]]) != 0){
          updated_value = Gibbs_update_omega(B=B, beta1[j], beta2[j], g=E[, j], temp, Zj=Z[i, list_map[[j]]], sigma2[i-1], sigma_g2[i-1, j], tau_g[j])
          idx = list_map[[j]][Z[i, list_map[[j]]] == 1]
          Omega[i, idx] = updated_value
        }
        else{
          Omega[i, list_map[[j]]] = 0
        }
        
        sigma_g2[i, j] = Gibbs_update_sigma2(E[, j], temp, Omega[i,list_map[[j]]], delta1, delta2)
        
        posterior[i] = posterior[i] + posterior_mechanistic_model(E[,j], M[ , list_map[[j]]], Omega[i,list_map[[j]] ], Zj=Z[i, list_map[[j]]], sigma_g2[i, j], tau_g[j], e, f, delta1, delta2)
        likelihood[i] = likelihood[i] + likelihood_mechanistic_model(E[,j], M[ , list_map[[j]]], Omega[i,list_map[[j]] ], sigma_g2[i, j])
      }
      
      ### update parameter in clinical model
      Ehat = methyl_block_product(M, list_map, Omega = Omega[i, ])
      X = cbind(E-Ehat, Ehat)
      
      for(k in 1:gammaTimes){
        gamma[(i-1)*gammaTimes + k , ] = MH_update_gamma(Y, X, gamma[(i-1)*gammaTimes + k-1, ], Hbar, tau, delta1, delta2, a, b, graph, rho=rho)
        
      }
      
      ### update beta
      err = Y - C %*% beta_c[i-1, ]
      updated_value = Gibbs_update_beta(err, X, gamma[i*gammaTimes, ], tau, sigma2[i-1])
      
      beta[i, gamma[i*gammaTimes, ] == 1] = updated_value
      beta[i, gamma[i*gammaTimes, ] != 1] = 0
      
      ### update beta_c
      err = Y - X %*% beta[i, ]
      beta_c[i, ] = Gibbs_update_beta_c(err, C, tau_c, sigma2[i-1] )
      
      ### update variance
      sigma2[i] = Gibbs_update_sigma2(Y, cbind(C, X), c(beta_c[i, ], beta[i, ]), delta1, delta2)
      
      posterior[i] = posterior[i] + posterior_clinic_model(Y, C, X, beta_c[i,], beta[i,], gamma[i*gammaTimes, ], sigma2[i], tau_c, tau, a, b, graph, delta1, delta2)
      likelihood[i] = likelihood[i] + likelihood_clinic_model(Y, C, X, beta_c[i,], beta[i,], sigma2[i])
      
      ### update missing data
      E_mis_bar = E[miss_met_sample, ] %*% beta[i, 1:p] 
      E_mis_hat = Ehat[miss_met_sample, ] %*% (beta[i, -(1:p)] - beta[i, 1:p])
      A_mis = Y[miss_met_sample] - C[miss_met_sample, ] %*% beta_c[i, ] - E_mis_hat - E_mis_bar
      for(j in 1:p)
      {
        beta_j1 = beta[i, 1:p][j]
        beta_j2 = beta[i, -(1:p)][j]
        A = A_mis + as.matrix(M[miss_met_sample, list_map[[j]] ]) %*% Omega[i, list_map[[j]] ]
        updated_value = methyl_update(A, E[miss_met_sample, j], as.matrix(M[miss_met_sample, list_map[[j]]]), Omega[i,list_map[[j]]], beta_j1, beta_j2, sigma2[i], sigma_g2[i, j])
        M[miss_met_sample, list_map[[j]] ] = updated_value
      }
    }
  }
  else if (n_exp_miss !=0 & n_met_miss ==0){
    # MCMC process
    for(i in 2:max_iters)
    {
      
      ### update parameter in mechanistic model
      B_tutle = Y - C %*% beta_c[i-1, ] - E %*% beta[i-1, 1:p]
      beta1 = beta[i-1, 1:p]
      beta2 = beta[i-1, -(1:p)]
      for(j in 1:p)
      {
        B = B_tutle - Ehat[, -j] %*% (beta2[-j] - beta1[-j])
        temp = as.matrix(M[, list_map[[j]]])
        Z[i, list_map[[j]]] = MH_update_Z(B=B, beta1[j], beta2[j], g=E[, j], temp, Zj=Z[i-1, list_map[[j]]], sigma2[i-1], sigma_g2[i-1, j], tau_g[j],e, f, rho=rho)
        
        if(sum(Z[i, list_map[[j]]]) != 0){
          updated_value = Gibbs_update_omega(B=B, beta1[j], beta2[j], g=E[, j], temp, Zj=Z[i, list_map[[j]]], sigma2[i-1], sigma_g2[i-1, j], tau_g[j])
          idx = list_map[[j]][Z[i, list_map[[j]]] == 1]
          Omega[i, idx] = updated_value
        }
        else{
          Omega[i, list_map[[j]]] = 0
        }
        
        sigma_g2[i, j] = Gibbs_update_sigma2(E[, j], temp, Omega[i,list_map[[j]]], delta1, delta2)
        
        posterior[i] = posterior[i] + posterior_mechanistic_model(E[,j], M[ , list_map[[j]]], Omega[i,list_map[[j]] ], Zj=Z[i, list_map[[j]]], sigma_g2[i, j], tau_g[j], e, f, delta1, delta2)
        likelihood[i] = likelihood[i] + likelihood_mechanistic_model(E[,j], M[ , list_map[[j]]], Omega[i,list_map[[j]] ], sigma_g2[i, j])
      }
      
      ### update parameter in clinical model
      Ehat = methyl_block_product(M, list_map, Omega = Omega[i, ])
      X = cbind(E-Ehat, Ehat)
      
      for(k in 1:gammaTimes){
        gamma[(i-1)*gammaTimes + k , ] = MH_update_gamma(Y, X, gamma[(i-1)*gammaTimes + k-1, ], Hbar, tau, delta1, delta2, a, b, graph, rho=rho)
        
      }
      
      ### update beta
      err = Y - C %*% beta_c[i-1, ]
      #print(sum(gamma[i*gammaTimes, ]))
      updated_value = Gibbs_update_beta(err, X, gamma[i*gammaTimes, ], tau, sigma2[i-1])
      
      beta[i, gamma[i*gammaTimes, ] == 1] = updated_value
      beta[i, gamma[i*gammaTimes, ] != 1] = 0
      
      ### update beta_c
      err = Y - X %*% beta[i, ]
      beta_c[i, ] = Gibbs_update_beta_c(err, C, tau_c, sigma2[i-1] )
      
      ### update variance
      sigma2[i] = Gibbs_update_sigma2(Y, cbind(C, X), c(beta_c[i, ], beta[i, ]), delta1, delta2)
      
      posterior[i] = posterior[i] + posterior_clinic_model(Y, C, X, beta_c[i,], beta[i,], gamma[i*gammaTimes, ], sigma2[i], tau_c, tau, a, b, graph, delta1, delta2)
      likelihood[i] = likelihood[i] + likelihood_clinic_model(Y, C, X, beta_c[i,], beta[i,], sigma2[i])
      
      ### update missing data
      E_mis_bar = E[miss_exp_sample, ] %*% beta[i, 1:p] 
      E_mis_hat = Ehat[miss_exp_sample, ] %*% (beta[i, -(1:p)] - beta[i, 1:p])
      A_mis = Y[miss_exp_sample] - C[miss_exp_sample, ] %*% beta_c[i, ] - E_mis_hat - E_mis_bar
      for(j in 1:p)
      {
        beta_j = beta[i,1:p][j]
        A = A_mis + E[miss_exp_sample,j] * beta_j
        
        E[miss_exp_sample, j] = mRNA_update(A=A, M=as.matrix(M[miss_exp_sample, list_map[[j]]]), beta_j, Omega[i, list_map[[j]] ], sigma2[i], sigma_g2[i,j])
      }
    }
  }
  else{
    # MCMC process
    for(i in 2:max_iters)
    {
      ### update parameter in mechanistic model
      B_tutle = Y - C %*% beta_c[i-1, ] - E %*% beta[i-1, 1:p]
      beta1 = beta[i-1, 1:p]
      beta2 = beta[i-1, -(1:p)]
      for(j in 1:p)
      {
        B = B_tutle - Ehat[, -j] %*% (beta2[-j] - beta1[-j])
        temp = as.matrix(M[, list_map[[j]]])
        Z[i, list_map[[j]]] = MH_update_Z(B=B, beta1[j], beta2[j], g=E[, j], temp, Zj=Z[i-1, list_map[[j]]], sigma2[i-1], sigma_g2[i-1, j], tau_g[j],e, f, rho=rho)
        
        if(sum(Z[i, list_map[[j]]]) != 0){
          updated_value = Gibbs_update_omega(B=B, beta1[j], beta2[j], g=E[, j], temp, Zj=Z[i, list_map[[j]]], sigma2[i-1], sigma_g2[i-1, j], tau_g[j])
          idx = list_map[[j]][Z[i, list_map[[j]]] == 1]
          Omega[i, idx] = updated_value
        }
        else{
          Omega[i, list_map[[j]]] = 0
        }
        
        sigma_g2[i, j] = Gibbs_update_sigma2(E[, j], temp, Omega[i,list_map[[j]]], delta1, delta2)
        
        posterior[i] = posterior[i] + posterior_mechanistic_model(E[,j], M[ , list_map[[j]]], Omega[i,list_map[[j]] ], Zj=Z[i, list_map[[j]]], sigma_g2[i, j], tau_g[j], e, f, delta1, delta2)
        likelihood[i] = likelihood[i] + likelihood_mechanistic_model(E[,j], M[ , list_map[[j]]], Omega[i,list_map[[j]] ], sigma_g2[i, j])
      }
      
      ### update parameter in clinical model
      Ehat = methyl_block_product(M, list_map, Omega = Omega[i, ])
      X = cbind(E-Ehat, Ehat)
      
      for(k in 1:gammaTimes){
        gamma[(i-1)*gammaTimes + k , ] = MH_update_gamma(Y, X, gamma[(i-1)*gammaTimes + k-1, ], Hbar, tau, delta1, delta2, a, b, graph, rho=rho)
        
      }
      
      ### update beta
      err = Y - C %*% beta_c[i-1, ]
      updated_value = Gibbs_update_beta(err, X, gamma[i*gammaTimes, ], tau, sigma2[i-1])
      
      beta[i, gamma[i*gammaTimes, ] == 1] = updated_value
      beta[i, gamma[i*gammaTimes, ] != 1] = 0
      
      ### update beta_c
      err = Y - X %*% beta[i, ]
      beta_c[i, ] = Gibbs_update_beta_c(err, C, tau_c, sigma2[i-1] )
      
      ### update variance
      sigma2[i] = Gibbs_update_sigma2(Y, cbind(C, X), c(beta_c[i, ], beta[i, ]), delta1, delta2)
      
      posterior[i] = posterior[i] + posterior_clinic_model(Y, C, X, beta_c[i,], beta[i,], gamma[i*gammaTimes, ], sigma2[i], tau_c, tau, a, b, graph, delta1, delta2)
      likelihood[i] = likelihood[i] + likelihood_clinic_model(Y, C, X, beta_c[i,], beta[i,], sigma2[i])
      
    }
  }
  
  return(list(beta_c = beta_c, 
              beta = beta, 
              gamma = gamma, 
              Omega = Omega, 
              Z = Z, 
              sigma2 = sigma2, 
              sigma_g2 = sigma_g2, 
              posterior=posterior, 
              likelihood = likelihood))
}