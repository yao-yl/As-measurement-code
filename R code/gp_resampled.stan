//  stan model for the resampled 271 wells
functions {
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order);
  vector build_b_spline(real[] t, real[] ext_knots, int ind, int order) {
    // INPUTS:
    //    t:          the points at which the b_spline is calculated
    //    ext_knots:  the set of extended knots
    //    ind:        the index of the b_spline
    //    order:      the order of the b-spline
    vector[size(t)] b_spline;
    vector[size(t)] w1 = rep_vector(0, size(t));
    vector[size(t)] w2 = rep_vector(0, size(t));
    if (order==1)
      for (i in 1:size(t)) // B-splines of order 1 are piece-wise constant
        b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind+1]);
    else {
      if (ext_knots[ind] != ext_knots[ind+order-1])
        w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t))) /
             (ext_knots[ind+order-1] - ext_knots[ind]);
      if (ext_knots[ind+1] != ext_knots[ind+order])
        w2 = 1 - (to_vector(t) - rep_vector(ext_knots[ind+1], size(t))) /
                 (ext_knots[ind+order] - ext_knots[ind+1]);
      // Calculating the B-spline recursively as linear interpolation of two lower-order splines
      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order-1) +
                 w2 .* build_b_spline(t, ext_knots, ind+1, order-1);
    }
    return b_spline;
  }
} 
data {
int N; // number of data
int num_knots;            // num of knots
vector[num_knots] knots;  // the sequence of knots
vector[N] y1;   // measurement in 00
vector[N] y2;    // measurement in 14
vector[N] y3;    // measurement in 15
row_vector[2] loc_X[N];
int spline_degree;        // the degree of spline (is equal to order - 1)
real rho;
real alpha;
real sigma;
int N_grid; // number of data
vector[N_grid] f_grid;
}
transformed data {
  real delta = 1e-9;
  int num_basis = num_knots + spline_degree - 1; // total number of B-splines
  matrix[N, N] L_K;
  matrix[N, N] K = cov_exp_quad(loc_X,alpha,rho) + diag_matrix(rep_vector(delta, N)); 
 matrix[num_basis, N_grid] B2; 
  L_K = cholesky_decompose(K);
  {
    vector[spline_degree + num_knots] ext_knots_temp;
    vector[2*spline_degree + num_knots] ext_knots; // set of extended knots
    ext_knots_temp = append_row(rep_vector(knots[1], spline_degree), knots);
    ext_knots = append_row(ext_knots_temp, rep_vector(knots[num_knots], spline_degree));
    for (ind in 1:num_basis)
       B2[ind,:] = to_row_vector(build_b_spline(to_array_1d(f_grid), to_array_1d(ext_knots), 
       ind, spline_degree + 1));
    B2[num_knots + spline_degree - 1, N_grid] = 1;
  }
}
parameters {
real a0; 
real mu;
real<lower=0> sigma_change;
row_vector[num_basis] a; 
vector[N] eta1;
vector[N] eta2;
}
transformed parameters{
vector[N] f1=  L_K * eta1+ mu;
vector[N] f2;
vector[N] change;
  {   
    matrix[num_basis, N] B; 
    vector[spline_degree + num_knots] ext_knots_temp;
    vector[2*spline_degree + num_knots] ext_knots; // set of extended knots
    ext_knots_temp = append_row(rep_vector(knots[1], spline_degree), knots);
    ext_knots = append_row(ext_knots_temp, rep_vector(knots[num_knots], spline_degree));
    for (ind in 1:num_basis)
       B[ind,:] = to_row_vector(build_b_spline(to_array_1d(f1), to_array_1d(ext_knots), 
       ind, spline_degree + 1));
    B[num_knots + spline_degree - 1, N] = 1;
    change= a0*f1 + to_vector(a*B); 
  }
  f2=f1+ change+ sigma_change* eta2 ;  
}
model{
  eta1~normal(0,1);
  eta2~normal(0,1);
  y1 ~ normal(f1, sigma);
  y2 ~ normal(f2, sigma);
  y3 ~ normal(f2, sigma);
	sigma_change~inv_gamma(3,3);
	a0 ~ normal(0, 0.5);
	a[1] ~ normal(0, 0.4);
	for (i in 1:(num_basis-1))
		a[i+1] ~ normal(a[i], 0.4);
	a~	normal(0,1.5);
	mu ~ normal(3.91, 0.5); // informed by sample mean;	
}

generated quantities{
  vector[N_grid] change_grid=a0*f_grid + to_vector(a*B2);
}

