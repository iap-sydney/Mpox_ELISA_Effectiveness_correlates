//
// This stan model fits the immunogenicity data
//

// The input data is a vector of reported GMTs (y_n), standard deviation (s_n),
// population (n_n) along with study index (study_ind), vaccine index (dose_ind)
// and formulation index (form_ind).
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

// Transform some data to have sum of squares (S_n) and population minus one (n_n_1)
transformed data{
  vector[N] S_n;
  vector[N] n_n_1;
  for (j in 1:N){
    n_n_1[j]=n_n[j]-1;
    S_n[j]=(n_n_1[j])*pow(s_n[j],2);
  }
}
// The parameters accepted in the model.
// mu_j=y_S
// ELISA effect is ignored here since it is not significant
parameters {
  vector[num_studies] mu_j;            // study effect
  real mu_f;                           // Formulation Effect
  vector[num_doses] mu_d;              // Mean for each vaccine schedule
  vector<lower=0>[num_doses] sigma_d;  // standard deviation for each schedule
  real<lower=0> sigma_s;               // Standard deviation of random effect
}

// Transform the parameters to define the sampling distributions
// sampling distributions are used as they are sufficient
transformed parameters{
  vector[N] mu_n;          // mean accounting for formulation, schedule and study
  vector[N] X_n;           // Chi-square statistic for sampling distribution
  vector[N] se_n;          // Standard error for sampling distribution
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

// Fitting the model to the data
model {
  y_n ~ normal(mu_n,se_n); // Fit GMT to sampling distribution 
    for (i in 1:N){
    target+=log(sigma_d[dose_ind[i]])*(-n_n[i]+1)-X_n[i]/2;  // Fit standard deviation
  }
  // Specify priors
  mu_d ~ normal (0,10);
  sigma_d ~ lognormal(0,10);
  mu_j ~ normal (0,sigma_s);
  sigma_s ~ cauchy(0,1);
  mu_f ~ normal(0,1);
}


