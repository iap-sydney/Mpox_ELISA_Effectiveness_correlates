//
// This Stan program contains the primary stan model for modelling the 
// decay of ELISA endpoint titers. 
//

// Define the decay functions used in our analysis
functions {
// Linear decay model used as a comparison
  real decay_linear(real time,real k1){
    return -k1*time;
  }
// 2 phase decay model implimented as our primary and best fitting model
  real decay_exp(real time,real k1,real k2, real f){
    return log10(f*exp(-k1*time)+(1-f)*exp(-k2*time));
  }
}
// The input data for our model.

data {
  int N;              // Number of Datapoints
  vector[N] y_n;      // Reported GMT for each datapoint
  vector[N] time;     // Time of GMT post peak (as defined in manuscript)
  vector[N] s_n;      // Standard deviation at each datapoint
  int n_n[N];         // Number of samples at each datapoint
  int num_doses;      // Number of vaccine regimens
  int num_form;       // Number of formulations used
  int num_studies;    // Number of studies included
  int study_ind[N];   // Study index for each datapoint
  int dose_ind[N];    // Index of dosing regimen for each datapoint
  vector[N] form_ind; // Index of formulations 
}

// Transform some of the data for various model inputs
transformed data{
  vector[N] S_n;     // Sum of squares at each data point
  vector[N] n_n_1;   // Number of participants minus 1
  for (j in 1:N){
    n_n_1[j]=n_n[j]-1;
    S_n[j]=(n_n_1[j])*pow(s_n[j],2);
  }
}
// The model parameters input into the model.
parameters {
  real<lower=0> decay1;                // Fast decay parameter
  real<lower=0,upper=decay1> decay2;   // Slow Decay parameter
  real<lower=0,upper=1> tc[num_doses]; // Proportion of fast decaying antibodies
  vector[num_studies] mu_j;            // Random effect of each study on GMT
  real mu_f;                           // Effect of formulation
  vector[num_doses] mu_d;              // GMT at time 0 for each regimen
  vector<lower=0>[num_doses] sigma_d;  // Standard deviation of titers per regimen
  real<lower=0> sigma_s;               // Standard deviation of random effect 
}

// Transform model parameters into those expressed in the model distributions
transformed parameters{
  vector[N] mu_n;                                 // Sampling Mean
  vector[N] X_n;                                  // Chi-Square statistic for ratio of sampling variance
  vector[N] se_n;                                 // Sampling Standard deviation
  for (i in 1:N){
    se_n[i]=sigma_d[dose_ind[i]]*pow(n_n[i],-0.5);
    X_n[i]=S_n[i]/pow(sigma_d[dose_ind[i]],2);
    if (form_ind[i]==2){
      mu_n[i]=mu_j[study_ind[i]]+mu_d[dose_ind[i]]+decay_exp(time[i],decay1,decay2,tc[dose_ind[i]]);
    } else {
    mu_n[i]=mu_j[study_ind[i]]+mu_d[dose_ind[i]]+mu_f+decay_exp(time[i],decay1,decay2,tc[dose_ind[i]]);
    } 
  }
}

// Input the Model parameters
model {
  y_n ~ normal(mu_n,se_n);     //Fit The observed GMTs to sampling distribution
  X_n ~ chi_square(n_n_1);     // Fit Observed standard deviations to sampling distribution
  // Define Priors
  mu_d ~ normal (0,10);        
  sigma_d ~ lognormal(0,10);
  mu_j ~ normal (0,sigma_s);
  sigma_s ~ cauchy(0,1);
  mu_f ~ normal(0,1);
  decay1~normal(0,1);
  decay2~normal(0,1);
}

// Generate R^2 fitted values for model comparison
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
