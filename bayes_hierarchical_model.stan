data {
  int<lower=0> N;//number of countries
  int<lower=0> CONT;//number of countinents
  int<lower=0> N1;//number of countries in NA
  int<lower=0> N2;//number of countries in EU
  int<lower=0> N3;//number of countries in AS
  int<lower=0> N4;//number of countries in AU
  int<lower=0> T;//time series length
  int<lower=0> C;//number of clusters
  int<lower=0> y[T, N, C]; //array of dimension T*N*C
  vector[C] p1dirich; // T=1 distribution for all countries
}


parameters {
  real<lower=-1,upper=1> alpha[C-1];
  vector<lower=0>[CONT] s0;
  cholesky_factor_corr[N1] L0_NA;
  cholesky_factor_corr[N2] L0_EU;
  cholesky_factor_corr[N3] L0_AS;
  vector[N1] eps_NA[T-1,C-1];
  vector[N2] eps_EU[T-1,C-1];
  vector[N3] eps_AS[T-1,C-1];
  real eps_AU[T-1,C-1];
  simplex[C] p1; // T=1 distribution for all countries
}

model {

  vector[N1] s0vec_NA;
  vector[N2] s0vec_EU;
  vector[N3] s0vec_AS;
  vector[N1] mu0vec_NA;
  vector[N2] mu0vec_EU;
  vector[N3] mu0vec_AS;
  real pexp[C-1];
  real lr_j[C-1];
  vector[C] p[T,N];
  vector[C] p_avg[T,N];
  
  
  p1 ~ dirichlet(p1dirich);
  s0[1] ~ cauchy(0, 0.5);
  s0[2] ~ cauchy(0, 0.5);
  s0[3] ~ cauchy(0, 0.5);
  s0[4] ~ cauchy(0, 0.5);
  
  for (i in 1:N){
    p[1,i,] = p1;
  }
  
  for (i in 1:N1){
    s0vec_NA[i] = s0[1];
    mu0vec_NA[i] = 0;
  }
  
  for (i in 1:N2){
    s0vec_EU[i] = s0[2];
    mu0vec_EU[i] = 0;
  }
  
  for (i in 1:N3){
    s0vec_AS[i] = s0[3];
    mu0vec_AS[i] = 0;
  }
  
  
  L0_NA ~ lkj_corr_cholesky(2.0);
  L0_EU ~ lkj_corr_cholesky(2.0);
  L0_AS ~ lkj_corr_cholesky(2.0);
  

  for (t in 2:T){
    
    for (j in 2:C) { 
      eps_NA[t-1,j-1] ~ multi_normal_cholesky(mu0vec_NA, diag_pre_multiply(s0vec_NA, L0_NA));
      eps_EU[t-1,j-1] ~ multi_normal_cholesky(mu0vec_EU, diag_pre_multiply(s0vec_EU, L0_EU));
      eps_AS[t-1,j-1] ~ multi_normal_cholesky(mu0vec_AS, diag_pre_multiply(s0vec_AS, L0_AS));
      eps_AU[t-1,j-1] ~ normal(0, s0[4]);
    }
    
    for (i in 1:N){
      
      if (i<=N1){
        for (j in 2:C) {
        lr_j[j-1] = log(p[t-1,i][j]/p[t-1,i][1]) + alpha[j-1] + eps_NA[t-1,j-1][i];
        }
        pexp = exp(lr_j);
      
        p[t,i][1] = 1/(sum(pexp)+1);
        for (j in 2:C) {
          p[t,i][j] = p[t,i][1] * pexp[j-1];
        }
      //print(p[t,i]);      
      }
      
      if (i>N1 && i<=N1+N2){
        for (j in 2:C) {
        lr_j[j-1] = log(p[t-1,i][j]/p[t-1,i][1]) + alpha[j-1] + eps_EU[t-1,j-1][i-N1];
        }
        pexp = exp(lr_j);
      
        p[t,i][1] = 1/(sum(pexp)+1);
        for (j in 2:C) {
          p[t,i][j] = p[t,i][1] * pexp[j-1];
        }
      //print(p[t,i]);      
      }
      
      if (i>N1+N2 && i<=N1+N2+N3){
        for (j in 2:C) {
        lr_j[j-1] = log(p[t-1,i][j]/p[t-1,i][1]) + alpha[j-1] + eps_AS[t-1,j-1][i-N1-N2];
        }
        pexp = exp(lr_j);
      
        p[t,i][1] = 1/(sum(pexp)+1);
        for (j in 2:C) {
          p[t,i][j] = p[t,i][1] * pexp[j-1];
        }
      //print(p[t,i]);      
      }
      
      if (i>N1+N2+N3){
        for (j in 2:C) {
        lr_j[j-1] = log(p[t-1,i][j]/p[t-1,i][1]) + alpha[j-1] + eps_AU[t-1,j-1];
        }
        pexp = exp(lr_j);
      
        p[t,i][1] = 1/(sum(pexp)+1);
        for (j in 2:C) {
          p[t,i][j] = p[t,i][1] * pexp[j-1];
        }
      //print(p[t,i]);     
      }
    
   }
      
  }
  
  

  for (i in 1:N){
    for (t in 2:T){
      
      //average p for first 14 days
      if (t < 15){
        for (j in 1:C){
          p_avg[t,i][j] = 0;
          for (k in 1:(t-1)){
            p_avg[t,i][j] = p_avg[t,i][j]+ p[t-k,i][j]/(t-1);
           //print(p_avg[t,i]);
          }
        }
      }
      
      //average p after first 14 days
      if (t >= 15){
        for (j in 1:C){
          p_avg[t,i][j] = 0;
          for (k in 1:14){
            p_avg[t,i][j] = p_avg[t,i][j]+ p[t-k,i][j]/14;
            //print(p_avg[t,i]);
          }
        }
      }
      
      y[t,i,1:C] ~ multinomial(p_avg[t,i]);
    }
  }
}


generated quantities {

  vector[N1] s0vec_NA;
  vector[N2] s0vec_EU;
  vector[N3] s0vec_AS;

  matrix[N1,N1] Omega_NA;
  matrix[N2,N2] Omega_EU;
  matrix[N3,N3] Omega_AS;
  matrix[N1,N1] Sigma_NA;
  matrix[N2,N2] Sigma_EU;
  matrix[N3,N3] Sigma_AS;

  for (i in 1:N1) {
    s0vec_NA[i] = s0[1];
  } 
  
  for (i in 1:N2) {
    s0vec_EU[i] = s0[2];
  }
  
  for (i in 1:N3) {
    s0vec_AS[i] = s0[3];
  }

  Omega_NA = multiply_lower_tri_self_transpose(L0_NA);
  Omega_EU = multiply_lower_tri_self_transpose(L0_EU);
  Omega_AS = multiply_lower_tri_self_transpose(L0_AS);
  Sigma_NA = quad_form_diag(Omega_NA, s0vec_NA);
  Sigma_EU = quad_form_diag(Omega_EU, s0vec_EU);
  Sigma_AS = quad_form_diag(Omega_AS, s0vec_AS);

}

