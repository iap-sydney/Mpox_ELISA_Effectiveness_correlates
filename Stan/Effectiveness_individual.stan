//
// This Stan program defines a model for comupting the effectiveness
// of multiple studies vaccines individually
//

// The input data consists of vectors of case and population. 
// Included is the disaggragation index (group_ind) and study index (study_ind)
data {
  int J;
  int num_cases[J];
  int Num_pop[J];
  int num_groups;
  int group_ind[J];
  int num_studies;
  int study_ind[J];
}

// This model consists of a base infection rate, r, for each group.
// 
parameters {
  real<lower=0,upper=1> r[num_groups];
  real A[num_studies];
}

transformed parameters{
  vector[num_studies] V;
  for (i in 1:num_studies){
    V[i]=log(inv_logit(A[i]));
  }
  vector<upper=1>[J] r_dose;
  for (i in 1:J){
    if (study_ind[i]<num_studies+1){
      r_dose[i]=inv_logit(logit(r[group_ind[i]])+V[study_ind[i]]);
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

