//
// This program fits a logistic model to the effectiveness  and immunogenicity
// data. A more in depth description of each section can be found in either 
// Effectiveness combined or Immunogenicity model for aspects relating to each
//

// Define the logistic function used to compute the expectation
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
  int num_vaccines;
  int Vaccine_ind[J];
  int num_vacc_studies;
  int vacc_study_ind[J];
  
  //miscellaneous sample
  int sample_size;
  vector[sample_size] sample;
}

//Transform the immunigenicity data
transformed data{
  vector[N] S_n;
  vector[N] n_n_1;
  for (j in 1:N){
    n_n_1[j]=n_n[j]-1;
    S_n[j]=(n_n_1[j])*pow(s_n[j],2);
  }
}
// The model parameters.
parameters {
  //immunogenicity model
  vector[num_studies] mu_j;
  real mu_f;
  vector[num_vaccines] mu_d;
  vector<lower=0>[num_vaccines] sigma_d;
  // Logistic Curve parameters
  real<lower=0> k;
  real A;
  // Effectiveness parameters
  real<lower=0> sigma_s;
  real<lower=0,upper=1> r[num_groups];
  real p_s[num_vacc_studies];
  real<lower=0> sigma;
}
// Transformed parameters
transformed parameters{
  //Immunogenicity Parameters
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
  // Protection using logistic function
  vector[num_vaccines] V;
  for (i in 1:num_vaccines){
    V[i]=log(1-eff(mu_d[i],k,A,sigma_d[i],sample));
  }
  // Effectiveness parameters
  vector<upper=1>[J] r_dose;
  for (i in 1:J){
    if (Vaccine_ind[i]<4){
      r_dose[i]=inv_logit(logit(r[group_index[i]])+V[Vaccine_ind[i]]+p_s[vacc_study_ind[i]]);
    }else{
      r_dose[i]=inv_logit(logit(r[group_index[i]]));
    }
  }
}

// Fitting the model.
model {
  // Immunogenicity data
  y_n ~ normal(mu_n,se_n);
  mu_d ~ normal (0,10);
  sigma_d ~ lognormal(0,10);
  mu_j ~ normal (0,sigma_s);
  sigma_s ~ cauchy(0,1);
  X_n ~ chi_square(n_n_1);
  mu_f ~ normal(0,1);
  // Effectiveness Data
  num_cases ~ binomial(Num_pop, r_dose);
  p_s ~ normal(0,sigma);
  sigma ~ cauchy(0,0.25);
  r ~ beta(1,1);
  // Priors on logistic curve
  k ~ normal(0,10);
  A ~ logistic(0,1);
}

