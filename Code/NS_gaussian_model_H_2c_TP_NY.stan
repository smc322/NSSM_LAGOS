
#Stan code for nonstationary model on lake total phosphorus in Northeastern region.
#Priors for other response variables and regions can be found in table 2 of supplementary material.


functions{
matrix kern_func(real th,real la){
matrix[2,2] G;
matrix[2,2] H;
G = [[la*cos(th), la*sin(th)], 
    [   -sin(th),    cos(th)]];
H = crossprod(G);
return H;
}}

data {
int<lower=1> N;       // Sample size
int<lower=1> S;       // Number of spatial locations (may be = N)
vector[N] y;          // Outcome variable
// Spatial coords
matrix[S,2] h_xy;     // vector of 2-D coordinates
// Mean function related data
int<lower=1> V_X;     // Number of mean function predictors 
matrix[N,V_X] X;      // Mean function design matrix
// variance-covariance design information
// int<lower=1> V_A;     // Number of spat-var alpha_sq predictors 
int<lower=1> V_L;     // Number of spat-var lambda predictors
int<lower=1> V_T;     // Number of spat-var theta predictors
// matrix[S,V_A] X_A;    // Spat-var alpha_sq design matrix       
matrix[S,V_L] X_L;    // Spat-var lambda design matrix
matrix[S,V_T] X_T;    // Spat-var theta design matrix
}

parameters {
// *** linear coefficient parameters ***
real Int_mu;  // Mean covariates parameters
vector[V_X] beta_mu;  // Mean covariates parameters
vector[V_L] beta_L;   // lambda covariates parameters
vector[V_T] beta_T;   // theta covariates parameters
// vector[V_A] beta_A;   // spatial var covariates parameters  

// *** Non-Centered hierarchical specifications *** 
vector[S] z_L;  //unscaled lambda hyper-parm z-score
vector[S] z_T;  //unscaled theta hyper-parm z-score
real<lower=0> theta_sig;    // hyper parm SDs of thetas and lambdas
real<lower=0> lambda_sig;

// *** sill, nugget, range ***
real<lower=0> tot_var;             //alpha^2 + nug^2
real<lower=0, upper=1> pct_nug;    //percent nugget
real<lower=0> rho;    //spatial range parameter 
}

transformed parameters{
matrix[S,S] Cov_NS; // Non-stationary corr matrix  
matrix[S,S] Q;      // Spat-var mahalanobis distance  
matrix[S,S] L_sp;   // Cholesky decomp matrix
real<lower=0> alpha_sq;
real<lower=0> nug_sq;

matrix[2,2] kern_mat_i;       //  2 x 2 aniso Kernel matrix  
matrix[2,2] kern_mat_j;       //  2 x 2 aniso Kernel matrix 
matrix[2,2] kern_avg_ij;      //  aniso Kernel matrix avg

vector[S] lambda = exp(X_L*beta_L + lambda_sig*z_L);  // non-centered Lambda
vector[S] theta_pct = inv_logit(X_T*beta_T + theta_sig*z_T); // non-centered Theta
// vector[S] theta_pct = Phi(X_T*beta_T + theta_sig*z_T); // non-centered Theta

// spatial var and nugget SD
alpha_sq = (1-pct_nug)*tot_var;
nug_sq =      pct_nug*tot_var;

// Off-diagonals of distance matrix Q & N-S Cov Mat
// For computational efficiency columns first
for (j in 1:(S-1)) {
  for (i in (j+1):S) {
// ith Kernel matrix ***
if(lambda[i]<=1.0)
  kern_mat_i = kern_func(pi()/2 + (pi()*theta_pct[i]), lambda[i]);
else
  kern_mat_i = kern_func(pi()*theta_pct[i], lambda[i]);
// jth Kernel matrix ***
if(lambda[j]<=1.0)
  kern_mat_j = kern_func(pi()/2 + (pi()*theta_pct[j]), lambda[j]);
else
  kern_mat_j = kern_func(pi()*theta_pct[j], lambda[j]);
// ij Kernel average ***
kern_avg_ij = 0.5*(kern_mat_i + kern_mat_j);

// Spat-varying Mahalanobis pseudo-distance matrix
Q[i,j] = quad_form(inverse(kern_avg_ij), to_vector(h_xy[i,]-h_xy[j,]));

// Use Q to compute the off-diag elements of Cov matrix
Cov_NS[i,j] = alpha_sq*(sqrt(sqrt(determinant(kern_mat_i)*determinant(kern_mat_j))) / 
                     sqrt(determinant(kern_avg_ij))) *
                     exp( -rho*sqrt(Q[i,j]) );  
// symmetry
Cov_NS[j,i] = Cov_NS[i,j]; 
}}
// diagonals
for(i in 1:S){
Cov_NS[i,i] = alpha_sq + nug_sq;}

// Take Cholesky decomp
L_sp = cholesky_decompose(Cov_NS); 
}

model {
vector[S] mu = Int_mu + X*beta_mu;     

target += normal_lpdf(Int_mu |2.8, 2.5);   // Intercept prior
target += normal_lpdf(beta_mu |0.0, 2.5);   // fixed effects prior
target += normal_lpdf(beta_L |0.0, 1.5);    // Lambda effects prior 
target += normal_lpdf(beta_T |0.0, 1.5);    // Theta effects prior 
// target += normal_lpdf(beta_A |0, 1);     // coefs for spat var SD  

target += normal_lpdf(z_L |0.0, 1.0);       // non-Centered unscaled Lambdas
target += normal_lpdf(z_T |0.0, 1.0);       // non-Centered unscaled Thetas
target += normal_lpdf(lambda_sig| 0.0, 0.2);// Lambdas' sigma parm
target += normal_lpdf(theta_sig| 0.0, 0.2); // Thetas' sigma parm

target += normal_lpdf(tot_var |0.0, 0.7);   // total var prior    
// target += beta_lpdf(pct_nug |2.0, 9.0);  // stronger nugget prior if needed
target += normal_lpdf(rho| 0.0, 1.8);       // weakly-informative prior for range

target += multi_normal_cholesky_lpdf(y | mu, L_sp);
}
generated quantities{
vector[S] mu_fit = Int_mu + X*beta_mu;
vector[S] Y_fit = multi_normal_cholesky_rng(mu_fit, L_sp);
real log_lik = multi_normal_cholesky_lpdf(y | mu_fit, L_sp);
real bayes_R2 = variance(Y_fit)/(variance(Y_fit) + variance(y - Y_fit));
}

