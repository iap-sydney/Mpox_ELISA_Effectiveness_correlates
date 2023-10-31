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
functions {
  real decay_linear(real time,real k1){
    return -1*k1*time;
  }
}
// The input data is a vector 'y' of length 'N'.
data {
  int N;
  vector[N] y_n;
  vector[N] time;
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
  real<lower=0> decay1;
  vector[num_studies] mu_j;
  real mu_f;
  vector[num_doses] mu_d;
  vector<lower=0>[num_doses] sigma_d;
  real<lower=0> var_s;
}

transformed parameters{
  vector[N] mu_n;
  vector[N] X_n;
  vector[N] se_n;
  for (i in 1:N){
    se_n[i]=sigma_d[dose_ind[i]]*pow(n_n[i],-0.5);
    X_n[i]=S_n[i]/pow(sigma_d[dose_ind[i]],2);
    if (form_ind[i]==1){
      mu_n[i]=mu_j[study_ind[i]]+mu_d[dose_ind[i]]+decay_linear(time[i],decay1);
    } else {
      mu_n[i]=mu_j[study_ind[i]]+mu_d[dose_ind[i]]+mu_f+decay_linear(time[i],decay1);
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
  mu_j ~ normal (0,pow(var_s,0.5));
  var_s ~ inv_gamma(1,1);
  X_n ~ chi_square(n_n_1);
  mu_f ~ normal(0,1);
  decay1~normal(0,1);
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N){
    log_lik[i]=normal_lpdf(y_n[i]|mu_n[i],se_n[i]);//+chi_square_lpdf(X_n[i]|n_n_1[i]);
  }
  real R;
  {
    vector[N] yres;
    yres=mu_n-y_n;
    R=variance(mu_n)/(variance(mu_n)+variance(yres));
  }
}
