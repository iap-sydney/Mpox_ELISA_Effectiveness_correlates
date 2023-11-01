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

// Input Data for the model
data {
  //Immunogenecity Data
  int N;               // Number of Immunogenicity datapoints
  vector[N] y_n;       // GMTs
  vector[N] s_n;       // Standard deviation
  int n_n[N];          // Number of samples in a datapoint
  int num_form;        // Numner of formulations
  int num_studies;     // Number of Immunogenicity studies
  int study_ind[N];    // Index of study
  int dose_ind[N];     // Index of vaccine
  vector[N] form_ind;  // Index of formulation
  
  //Effectiveness case data
  int J;                 // Number of Efficacy datapoints
  int num_cases[J];      // Number of cases 
  int Num_pop[J];        // Population size
  int num_groups;        // Number of distinct groups in efficacy
  int group_index[J];    // Group Index
  int num_vaccines;      // Number of Vaccines
  int Vaccine_ind[J];    // Index of vaccine used
  int num_vacc_studies;  // Number of Efficacy Studies
  int vacc_study_ind[J]; // Index of efficacy study
  
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
  // Priors on Immunogenicity
  sigma_d ~ lognormal(0,10);
  mu_j ~ normal (0,sigma_s);
  sigma_s ~ cauchy(0,1);
  X_n ~ chi_square(n_n_1);
  mu_f ~ normal(0,1);
  // Effectiveness Data
  num_cases ~ binomial(Num_pop, r_dose);
  // Prior on Efficacy data
  p_s ~ normal(0,sigma);
  sigma ~ cauchy(0,0.25);
  r ~ beta(1,1);
  // Priors on logistic curve
  k ~ normal(0,10);
  A ~ logistic(0,1);
}

