//
// This Stan program defines a model for comupting the combined effectiveness
// accross different studies
//

// The input data consists of vectors of case and population. 
// Included is the disaggragation index (group_ind),
// vaccine inex (vaccine_ind) and study index (study_ind)
data {
  int J;
  int num_cases[J];
  int Num_pop[J];
  int num_groups;
  int group_index[J];
  int num_vaccines;
  int Vaccine_ind[J];
  int num_vacc_studies;
  int vacc_study_ind[J];
}

// This model consists of a base infection rate, r, for each group.
// And an auxillary parameter,A, for the effectiveness fo consistency in analysis
// A=exp(R_v)
parameters {
  real<lower=0,upper=1> r[num_groups];
  real p_s[num_vacc_studies];
  real<lower=0,upper=1> A[num_vaccines];
  real<lower=0> sigma;
}

//Transform A parameters into model parameters
//V=R_v
//r_dose=r_i,v,s
transformed parameters{
  vector[num_vaccines] V;
  for (i in 1:num_vaccines){
    V[i]=log(A[i]);
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

// Fit the binomial model with relevant priors.
model {
  num_cases ~ binomial(Num_pop, r_dose);
  p_s ~ normal(0,sigma);
  A ~ uniform(0,1);
  sigma ~ cauchy(0,0.25);
  r ~ beta(1,1);
}



