data{
  int<lower=0> n1;  // the number of wells in survey 1,
  int<lower=0> n2; // the number of wells in survey 2,
  int w; //  the number of discrete levels,
  int num_basis; //  the number of knots,
  matrix[n1, num_basis] B1;   // 2-D basis functions at x1
  matrix[n2, num_basis] B2;   // 2-D basis functions at x2
  matrix[n2, num_basis] Bd2;   // Laplacian of 2-D basis functions at x2
  int<lower=0> B1_m;  // the number non-zero component
  int<lower=0> B2_m;  // the number non-zero component
  int<lower=0> Bd2_m;  // the number non-zero component
  real  individual_y1[n1];
  vector[n1]  d1;  //depth  at the first survey  
  vector[n2]  d2; //depth  at the second survey  
  int <lower=0, upper=w> z2[n2]; // second survey response (categorical)
  ordered[w-1] c; // pre-learned logistic regression 
  real logistic_beta; // pre-learned logistic regression 
}

transformed data{
  vector[B1_m] w1=csr_extract_w(B1);
  vector[B2_m] w2=csr_extract_w(B2);
  vector[Bd2_m] wd2=csr_extract_w(Bd2);
  int v1[B1_m];
  int v2[B2_m];
  int vd2[Bd2_m];
  int u1[n1+1];
  int u2[n2+1];
  int ud2[n2+1];
  v1=csr_extract_v(B1);
  v2=csr_extract_v(B2);
  vd2=csr_extract_v(Bd2);
  u1=csr_extract_u(B1);
  u2=csr_extract_u(B2);
  ud2=csr_extract_u(Bd2);
}

parameters{
  real<multiplier=0.3> beta_l;
  real<multiplier=0.3> beta_y;
  real<multiplier=0.5> beta_y_log;
  real<lower=0> sigma;
  real<lower=0> sigma2;
  real beta_itercept;
  vector[num_basis] a; 
  real a0;
  vector [n2]   obs_error_2;
  real beta_depth;
}
transformed parameters{
  vector[n2]  laplician=csr_matrix_times_vector(n2, num_basis, wd2, vd2, ud2, a );     //  Bd2*a   ;
  vector[n2]  eta2;
  vector[n2]  y2_individual;
  vector[n2]  y21=csr_matrix_times_vector(n2, num_basis, w2, v2, u2, a )+a0 + beta_depth*d2;       // B2*a + a0; 
  vector[n1]  y1_hat= csr_matrix_times_vector(n1, num_basis, w1, v1, u1, a )+a0+ beta_depth*d1;       // B1*a + a0;
  {
      vector[n2]  y21_exp= exp(y21/2);
      y2_individual = beta_itercept + y21 + (beta_y*y21_exp+ beta_y_log*y21)  .* laplician+ beta_l* laplician + (sigma+sigma2)* obs_error_2; 
  }
  eta2=logistic_beta*y2_individual;
}
model{
  individual_y1 ~ normal(y1_hat, sigma); 
  z2~ordered_logistic(eta2, c);
  obs_error_2 ~ normal (0,1);
  a~normal (0,0.5);
  sigma~inv_gamma(5,5);
  sigma2~inv_gamma(3,3);
  a0~normal (4,2);
  beta_y~normal (0,0.2);
  beta_l~normal (0,0.1);
  beta_y_log~normal (0,0.5);
  beta_itercept~normal (0,1);
  beta_depth~normal(0,0.6);
}

