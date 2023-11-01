//
// This Stan program defines an extended decay model with different decay rates 
// for each regimen used. Where code overlaps with Decay_model.stan, description
// is provided in Decat_model.stan
//

functions {
  real decay_fun(real time,real k1,real k2, real f){
    return log10(f*exp(-k1*time)+(1-f)*exp(-k2*time));
  }
}

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

parameters {
  real<lower=0> decay1[num_doses];              // Decay rate per vaccination schedule
  real<lower=0,upper=decay1> decay2[num_doses];
  real<lower=0,upper=1> tc[num_doses];
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
      mu_n[i]=mu_j[study_ind[i]]+mu_d[dose_ind[i]]+decay_fun(time[i],decay1[dose_ind[i]],decay2[dose_ind[i]],tc[dose_ind[i]]);
    } else {
      mu_n[i]=mu_j[study_ind[i]]+mu_d[dose_ind[i]]+mu_f+decay_fun(time[i],decay1[dose_ind[i]],decay2[dose_ind[i]],tc[dose_ind[i]]);
    } 
  }
}


model {
  y_n ~ normal(mu_n,se_n);
  mu_d ~ normal (0,10);
  sigma_d ~ lognormal(0,10);
  mu_j ~ normal (0,sigma_s);
  sigma_s ~ cauchy(0,1);
  X_n ~ chi_square(n_n_1);
  mu_f ~ normal(0,1);
  decay1~normal(0,1);
  decay2~normal(0,0.1);
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N){
    log_lik[i]=normal_lpdf(y_n[i]|mu_n[i],se_n[i]);// + chi_square_lpdf(X_n[i]|n_n_1[i]);
  }
    real R;
  {
    vector[N] yres;
    yres=mu_n-y_n;
    R=variance(mu_n)/(variance(mu_n)+variance(yres));
  }
}
