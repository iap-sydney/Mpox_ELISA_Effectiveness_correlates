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
// And an auxillary parameter,A, for the effectiveness fo consistency in analysis
// A=exp(R_v)
parameters {
  real<lower=0,upper=1> r[num_groups];
  real<lower=0,upper=1> A[num_studies];
}

//Transform A parameters into model parameters
//V=R_v
//r_dose=r_i,v,s
transformed parameters{
  vector[num_studies] V;
  for (i in 1:num_studies){
    V[i]=log(A[i]);
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

// Fit the binomial model with relevant priors.
model {
  num_cases ~ binomial(Num_pop, r_dose);
  A ~ uniform(0,1);
  r ~ beta(1,1);
}
