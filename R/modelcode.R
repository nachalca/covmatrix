sim.iw = "
data {
  int <lower=0> N;
  int <lower=0> k;  
  matrix[k,k] R;
  vector[k] y[N];
  vector[k] mu0;
}
parameters {
  //  vector[k] mu;
  cov_matrix[k] Sigma;
}
transformed parameters {
  real s1;
  real s2;
  real rho;
  s1 <- sqrt(Sigma[1,1]);
  s2 <- sqrt(Sigma[2,2]);
  rho <- Sigma[1,2]/(s1*s2);
}
model {
//  for ( i in 1:k)  mu[i] ~ normal(0, 100);
  Sigma ~ inv_wishart(k+1, R);
  for (n in 1:N) y[n] ~ multi_normal(mu0, Sigma);
}
"

#---------------------
# Scaled data IW
# TO DO
#---------------------------

# SIW model
# the actual covariance matrix in siw model would be Sigma=D*Q*D
sim.siw = "
data {
  int <lower=0> N;
  int <lower=0> k;
  matrix[k,k] R;
  vector[k] y[N];
  vector[k] mu0;
}
parameters {
//  vector[k] mu;
  cov_matrix[k] Q;
  vector[k] xi;
}
transformed parameters {
  matrix[k,k] D;
  vector[k] delta;

  real s1;
  real s2;
  real rho;
  for (i in 1:k) delta[i] <- exp( xi[i] );

D <- diag_matrix(delta);
s1 <- sqrt( Q[1,1]*delta[1]*delta[1]);
s2 <- sqrt(Q[2,2]*delta[2]*delta[2]);
rho <-  (Q[1,2]*delta[2]*delta[1]) /(s1*s2);
}
model {
  matrix[k,k] L;
  matrix[k,k] A;
  Q ~ inv_wishart(k+1, R);
  L <- cholesky_decompose(Q);
  A <- D*L;
    for ( i in 1:k) {
    //    mu[i] ~ normal(0, 100);
        xi[i] ~ normal(0, 100);
    }
  for (n in 1:N) y[n] ~ multi_normal_cholesky(mu0, A);
}
"
#------------------------

sim.ssig = "
data {
  int <lower=0> N;
  int <lower=0> k;
  matrix[k,k] R;
  vector[k] y[N];
  vector[k] mu0;
}
parameters {
  //  vector[k] mu;
  cov_matrix[k] Q1;
  vector[k] xi;
}
transformed parameters {
  matrix[k,k] L;
  corr_matrix[k] Rho; 
  cov_matrix[k] Sigma;
  vector<lower=0>[k] delta;
  vector<lower=0>[k] delta1;
  real s1;
  real s2;
  real rho;
  // Rho is the correlation matrix prior, start with a Q1 ~ IW() and its transformed into
  // a correlation matrix with D1*Q1*D1, wehre D1<-diag(delta1), is done with for loops
  
  for (i in 1:k) delta1[i] <- 1/sqrt(Q1[i,i]);
  
  for (n in 1:k) {
    for (m in 1:n) {
      Rho[m,n] <- delta1[m] * delta1[n] * Q1[m,n]; 
    }
  }
  
  for (n in 1:k) {
    for (m in (n+1):k) {
      Rho[m,n] <- Rho[n,m];
    }
  } 
  
  // compute covariance matrix as: Sigma = D*Q*D, where D = diag(delta) 
  //for (i in 1:k)  delta[i] <- exp( xi[i] );
  for (i in 1:k)  delta[i] <-  sqrt(xi[i]) ;
  for (n in 1:k) {
    for (m in 1:n) {
      Sigma[m,n] <- delta[m] * delta[n] * Rho[m,n]; 
    }
  }
  for (n in 1:k) {
    for (m in (n+1):k) {
      Sigma[m,n] <- Sigma[n,m];
    }
  }
  s1 <- sqrt( Sigma[1,1]);
  s2 <- sqrt(Sigma[2,2]) ;
  rho <-  Sigma[1,2] /(s1*s2);
}
model {
  Q1 ~ inv_wishart(k+1, R);
  for ( i in 1:k) {
    //      mu[i] ~ normal(0, 100);
     // xi[i] ~ normal(0, log(100));
  xi[i] ~ inv_gamma(1, .5); // implied variance prior from IW(d+1, I)
  }
  for (n in 1:N) y[n] ~ multi_normal(mu0, Sigma);
}
"

