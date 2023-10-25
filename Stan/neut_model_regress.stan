//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int N;
  vector[N] y_n;
  vector[N] s_n;
  int n_n[N];
  int num_doses;
  int num_form;
  int num_studies;
  int study_ind[N];
  int dose_ind[N];
  vector[N] form_ind;
}
transformed data{
  vector[N] S_n;
  vector[N] n_n_1;
  for (j in 1:N){
    n_n_1[j]=n_n[j]-1;
    S_n[j]=(n_n_1[j])*pow(s_n[j],2);
  }
}
// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[num_studies] mu_j;
  real mu_f;
  vector[num_doses] mu_d;
  vector<lower=0>[num_doses] sigma_d;
  real<lower=0> sigma_s;
}

transformed parameters{
  vector[N] mu_n;
  vector[N] X_n;
  vector[N] se_n;
  for (i in 1:N){
    se_n[i]=sigma_d[dose_ind[i]]*pow(n_n[i],-0.5);
    X_n[i]=S_n[i]/pow(sigma_d[dose_ind[i]],2);
    if (form_ind[i]==1){
      mu_n[i]=mu_j[study_ind[i]]+mu_d[dose_ind[i]]+mu_f;
    } else {
      mu_n[i]=mu_j[study_ind[i]]+mu_d[dose_ind[i]];
    }
  }
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  y_n ~ normal(mu_n,se_n);
  mu_d ~ normal (0,10);
  sigma_d ~ lognormal(0,10);
  mu_j ~ normal (0,sigma_s);
  sigma_s ~ cauchy(0,1);
  //X_n ~ chi_square(n_n_1);
  for (i in 1:N){
    target+=log(sigma_d[dose_ind[i]])*(-n_n[i]+1)-X_n[i]/2;
  }
  mu_f ~ normal(0,1);
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N){
    log_lik[i]=normal_lpdf(y_n[i]|mu_n[i],se_n[i])+chi_square_lpdf(X_n[i]|n_n_1[i]);
  }
  real Tx;
  real T_x;
  real T_xrep;
  {
    vector[N] sim_means;
    vector[N] test_mu;
      for(i in 1:N){
        if (form_ind[i]==1){
          test_mu[i]=mu_d[dose_ind[i]]+mu_f;
        } else {
          test_mu[i]=mu_d[dose_ind[i]];
        }
        sim_means[i]=normal_rng(test_mu[i],se_n[i]);
      }
    Tx=sum(pow(sim_means,2));
    for (i in 1:N){
      sim_means[i]=normal_rng(mu_n[i],se_n[i]);
    }
    
    
    T_xrep=sum(pow(sim_means-mu_n,2));
    T_x=sum(pow(y_n-mu_n,2));
  }
}

