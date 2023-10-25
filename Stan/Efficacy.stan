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
  int J;
  int num_cases[J];
  int Num_pop[J];
  int num_groups;
  int group_ind[J];
  int num_vaccines;
  int Vaccine_ind[J];
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real<lower=0,upper=1> r[num_groups];
  real<upper=0> A[num_vaccines];
}

transformed parameters{
  vector[num_vaccines] V;
  for (i in 1:num_vaccines){
    V[i]=log(inv_logit(A[i]));
  }
  vector<upper=1>[J] r_dose;
  for (i in 1:J){
    if (Vaccine_ind[i]<num_vaccines+1){
      r_dose[i]=inv_logit(logit(r[group_ind[i]])+V[Vaccine_ind[i]]);
    } else {
      r_dose[i]=r[group_ind[i]];
    }
  }
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  num_cases ~ binomial(Num_pop, r_dose);
  A ~ logistic(0,1);
  r ~ beta(1,1);
}

generated quantities {
  vector[J] log_lik;
  {
    for (i in 1:J){
      log_lik[i]=binomial_lpmf(num_cases[i]|Num_pop[i],r_dose[i]);
    }
  }
}

