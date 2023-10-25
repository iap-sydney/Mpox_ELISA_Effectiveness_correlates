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
functions{
  real eff(real mu,real k,real A,real sigma,vector sample){
        return mean(pow(1+exp((sample*sigma+mu-2.5)*(-k)-A),-1));
  }
}

// The input data is a vector 'y' of length 'N'.
data {
  //Immunogenecity Data
  int N; #Number of studies
  vector[N] y_n;
  vector[N] s_n;
  int n_n[N];
  int num_form;
  int num_studies;
  int study_ind[N];
  int dose_ind[N];
  vector[N] form_ind;
  
  //Effectiveness case data
  int J;
  int num_cases[J];
  int Num_pop[J];
  int num_groups;
  int group_index[J];
  int num_2dose;
  int num_vaccines;
  int Vaccine_ind[J];
  int num_vacc_studies;
  int vacc_study_ind[J];
  vector[5000] sample;
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
  vector[num_vaccines] mu_d;
  vector<lower=0>[num_vaccines] sigma_d;
  real<lower=0> sigma_s;
  real<lower=0> k;
  real A;
  real<lower=0,upper=1> r[num_groups];
  real p_s[num_vacc_studies];
  real<lower=0> sigma;
}

transformed parameters{
  vector[N] mu_n;
  vector[N] X_n;
  vector[N] se_n;
  for (i in 1:N){
    se_n[i]=sigma_d[dose_ind[i]]*pow(n_n[i],-0.5);
    X_n[i]=S_n[i]/pow(sigma_d[dose_ind[i]],2);
    if (form_ind[i]==3){
      mu_n[i]=mu_j[study_ind[i]]+mu_d[dose_ind[i]];
    } else if (form_ind[i]==2){
      mu_n[i]=mu_j[study_ind[i]]+mu_d[dose_ind[i]]+mu_f;
    } else {
      mu_n[i]=mu_j[study_ind[i]]+mu_d[dose_ind[i]];
    }
  }
  vector[num_vaccines] V;
  //V=log(1-pow(1+exp(-(mu_d-2.5)*k-A),-1));
  for (i in 1:num_vaccines){
    V[i]=log(1-eff(mu_d[i],k,A,sigma_d[i],sample));
  }
  vector<upper=1>[J] r_dose;
  for (i in 1:J){
    if (Vaccine_ind[i]<4){
      r_dose[i]=inv_logit(logit(r[group_index[i]])+V[Vaccine_ind[i]]+p_s[vacc_study_ind[i]]);
    }else{
      r_dose[i]=inv_logit(logit(r[group_index[i]]));
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
  X_n ~ chi_square(n_n_1);
  mu_f ~ normal(0,1);
  num_cases ~ binomial(Num_pop, r_dose);
  p_s ~ normal(0,sigma);
  k ~ normal(0,10);
  A ~ logistic(0,1);
  sigma ~ cauchy(0,0.25);
  r ~ beta(1,1);
}


generated quantities {
  vector[J] log_lik;
  //vector[num_vacc_studies] log_lik;
  {
    for (i in 1:J){
      log_lik[i]=binomial_lpmf(num_cases[i]|Num_pop[i],r_dose[i]);
    }
    //for (i in 1:num_vacc_studies){
    //log_lik[i]=normal_lpdf(p_s[i]|0,sigma);
    //}
  }
  int sim_cases[J];
  sim_cases=binomial_rng(Num_pop,r_dose);
}
