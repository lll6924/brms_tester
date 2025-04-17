// generated with brms 2.22.2
functions {
 /* compute correlated group-level effects
  * Args:
  *   z: matrix of unscaled group-level effects
  *   SD: vector of standard deviation parameters
  *   L: cholesky factor correlation matrix
  * Returns:
  *   matrix of scaled group-level effects
  */
  matrix scale_r_cor(matrix z, vector SD, matrix L) {
    // r is stored in another dimension order than z
    return transpose(diag_pre_multiply(SD, L) * z);
  }


   /* collect elements based on indices
   * Args:
   *   elements: response value
   *   J: the index of each element
   * Returns:
   *   a vector of size max(J)
   */
    vector bincount(vector elements, array[] int J) {
        vector[max(J)] res = rep_vector(0, max(J));
        for (n in 1:size(J)) {
          res[J[n]] = res[J[n]] + elements[n];
        }
        return res;
    }

   /* batched inversion of spd matrices
   * Args:
   *   ar: batched spd matrices as a 3d array
   * Returns:
   *   batched spd matrices (inversions of input) as a 3d array
   */
    array[,,] real batch_inverse_spd(array[,,] real ar){
        array[size(ar),size(ar[1]),size(ar[1])] real res;
        for (n in 1:size(ar)) {
            res[n] = to_array_2d(inverse_spd(to_matrix(ar[n])));
        }
        return res;
    }

   /* sum of log determinants of batched spd matrices
   * Args:
   *   ar: batched spd matrices as a 3d array
   * Returns:
   *   a real number indicating the summation of log determinants
   */
    real batch_log_determinant_spd(array[,,] real ar){
        real res = 0;
        for (n in 1:size(ar)) {
            res = res + log_determinant_spd(to_matrix(ar[n]));
        }
        return res;
    }

   /* log pdf of a normal_id_glm mixed-effect model with varying noise scales where one term with a single data vector is marginalized out
   * Args:
   *   Y: a vector of the response variable
   *   Xc: the data matrix in the linear regression part
   *   mu: the mean of the mixed-effect parts that are not marginalized
   *   b: the coefficients in the linear regression part
   *   sigma: a vector of the observation noise scales
   *   J: an array of the grouping indices of the marginalized effect
   *   N: number of groups in the marginalized effect
   *   M: number of data vectors in the marginalized effect, M=1 in this function
   *   sd: a vector of the noise scales for the marginalized effect
   *   Z: the data vector
   * Returns:
   *   a real number indicating the log pdf
   */
    real normal_id_glm_marginalized_lpdf(vector Y, matrix Xc, vector mu, vector b, vector sigma, array[] int J, int N, int M, vector sd, vector Z){
        int n_obs = size(Y);
        vector[n_obs] z = Y - mu - Xc * b;
        real s_u = sd[1] * sd[1];
        vector[n_obs] s_y = sigma .* sigma;
        vector[N] F = 1 / s_u + bincount(Z .* Z./ s_y, J) ;
        vector[N] x = bincount(z.* Z ./s_y, J);
        real a = sum(log(F)) + log(s_u) * N + sum(log(s_y));
        real c = sum(z .* z ./ s_y) - sum(x .* x ./ F);
        return - (a + c) / 2;
        //matrix[1,1] L = rep_matrix(1,1,1);
        //array[1] vector[size(Y)] ZZ;
        //ZZ[1] = Z;
        //return multi_normal_id_glm_marginalized_lpdf(Y| Xc, mu, b, sigma, J, N, M, sd, L, ZZ);
    }

   /* log pdf of a normal_id_glm mixed-effect model with fixed noise scales where one term with a single data vector is marginalized out
   * Args:
   *   Y: a vector of the response variable
   *   Xc: the data matrix in the linear regression part
   *   mu: the mean of the mixed-effect parts that are not marginalized
   *   b: the coefficients in the linear regression part
   *   sigma: the fixed observation noise scale
   *   J: an array of the grouping indices of the marginalized effect
   *   N: number of groups in the marginalized effect
   *   M: number of data vectors in the marginalized effect, M=1 in this function
   *   sd: a vector of the noise scales for the marginalized effect
   *   Z: the data vector
   * Returns:
   *   a real number indicating the log pdf
   */
    real normal_id_glm_marginalized_lpdf(vector Y, matrix Xc, vector mu, vector b, real sigma, array[] int J, int N, int M, vector sd, vector Z){
        int n_obs = size(Y);
        vector[n_obs] sigma2 = rep_vector(sigma, n_obs);
        return normal_id_glm_marginalized_lpdf(Y| Xc, mu, b, sigma2, J, N, M, sd, Z);
    }

   /* log pdf of a normal_id_glm mixed-effect model with varying noise scales where one term with multiple data vectors is marginalized out
   * Args:
   *   Y: a vector of the response variable
   *   Xc: the data matrix in the linear regression part
   *   mu: the mean of the mixed-effect parts that are not marginalized
   *   b: the coefficients in the linear regression part
   *   sigma: a vector of the observation noise scales
   *   J: an array of the grouping indices of the marginalized effect
   *   N: number of groups in the marginalized effect
   *   M: number of data vectors in the marginalized effect, M>1 in this function
   *   sd: a vector of the noise scales for the marginalized effect
   *   L: Cholesky factor of the correlation matrix for the marginalized effect
   *   Z: an array of the data vectors
   * Returns:
   *   a real number indicating the log pdf
   */
    real normal_id_glm_multi_marginalized_lpdf(vector Y, matrix Xc, vector mu, vector b, vector sigma, array[] int J, int N, int M, vector sd, matrix L, array[] vector Z){
        int n_obs = size(Y);
        vector[n_obs] z = Y - mu - Xc * b;
        matrix[M,M] s_u = multiply_lower_tri_self_transpose(diag_pre_multiply(sd, L));
        matrix[M,M] inv_s_u = inverse_spd(s_u);
        vector[n_obs] s_y = sigma .* sigma;
        array[N,M,M] real F;
        for (m1 in 1:M){
            for (m2 in 1:m1){
                F[:,m1,m2] = to_array_1d(bincount(Z[m1] .* Z[m2]./ s_y, J) );
                F[:,m2,m1] = F[:,m1,m2];
            }
        }
        for (n in 1:N){
            F[n] = to_array_2d(to_matrix(F[n]) + inv_s_u);
        }
        array[N,M,M] real F_inv = batch_inverse_spd(F);
        array[N,M] real x;
        for (m in 1:M){
            x[:,m] = to_array_1d(bincount(z .* Z[m] ./ s_y, J));
        }
        real a = batch_log_determinant_spd(F) + log_determinant_spd(s_u) * N + sum(log(s_y));
        real c = sum(z .* z ./ s_y);
        for (n in 1:N){
            c = c - quad_form(to_matrix(F_inv[n]), to_vector(x[n]));
        }
        return - (a + c) / 2;
    }

   /* log pdf of a normal_id_glm mixed-effect model with fixed noise scales where one term with multiple data vectors is marginalized out
   * Args:
   *   Y: a vector of the response variable
   *   Xc: the data matrix in the linear regression part
   *   mu: the mean of the mixed-effect parts that are not marginalized
   *   b: the coefficients in the linear regression part
   *   sigma: the fixed observation noise scale
   *   J: an array of the grouping indices of the marginalized effect
   *   N: number of groups in the marginalized effect
   *   M: number of data vectors in the marginalized effect, M>1 in this function
   *   sd: a vector of the noise scales for the marginalized effect
   *   L: Cholesky factor of the correlation matrix for the marginalized effect
   *   Z: an array of the data vectors
   * Returns:
   *   a real number indicating the log pdf
   */
    real normal_id_glm_multi_marginalized_lpdf(vector Y, matrix Xc, vector mu, vector b, real sigma, array[] int J, int N, int M, vector sd, matrix L, array[] vector Z){
        int n_obs = size(Y);
        vector[n_obs] sigma2 = rep_vector(sigma, n_obs);
        return normal_id_glm_multi_marginalized_lpdf(Y| Xc, mu, b, sigma2, J, N, M, sd, L, Z);
    }

   /* recovery of the marginalized effect in a normal_id_glm mixed-effect model with varying noise scales where one term with a single data vector is marginalized out
   * Args:
   *   Y: a vector of the response variable
   *   Xc: the data matrix in the linear regression part
   *   mu: the mean of the mixed-effect parts that are not marginalized
   *   b: the coefficients in the linear regression part
   *   sigma: a vector of the observation noise scales
   *   J: an array of the grouping indices of the marginalized effect
   *   N: number of groups in the marginalized effect
   *   M: number of data vectors in the marginalized effect, M=1 in this function
   *   sd: a vector of the noise scales for the marginalized effect
   *   Z: the data vector
   * Returns:
   *   an array of vectors indicating the marginalized effects
   */
    array[] vector normal_id_glm_marginalized_recover_rng(vector Y, matrix Xc, vector mu, vector b, vector sigma, array[] int J, int N, int M, vector sd, vector Z){
        int n_obs = size(Y);
        vector[n_obs] z = Y - mu - Xc * b;
        real s_u = sd[1] * sd[1];
        vector[n_obs] s_y = sigma .* sigma;
        vector[N] G = bincount(Z .* Z./ s_y, J);
        vector[N] F = 1 / s_u + G;
        vector[N] x = bincount(z .* Z ./s_y, J);
        vector[N] Mu = s_u * (1 - G ./ F) .* x;
        vector[N] L = sqrt(1/F);
        array[M] vector[N] ans;
        for (n in 1:N) {
            ans[1,n] = normal_rng(Mu[n], L[n]);
        }
        return ans;
    }

   /* recovery of the marginalized effect in a normal_id_glm mixed-effect model with fixed noise scales where one term with a single data vector is marginalized out
   * Args:
   *   Y: a vector of the response variable
   *   Xc: the data matrix in the linear regression part
   *   mu: the mean of the mixed-effect parts that are not marginalized
   *   b: the coefficients in the linear regression part
   *   sigma: the fixed observation noise scales
   *   J: an array of the grouping indices of the marginalized effect
   *   N: number of groups in the marginalized effect
   *   M: number of data vectors in the marginalized effect, M=1 in this function
   *   sd: a vector of the noise scales for the marginalized effect
   *   Z: the data vector
   * Returns:
   *   an array of vectors indicating the marginalized effects
   */
    array[] vector normal_id_glm_marginalized_recover_rng(vector Y, matrix Xc, vector mu, vector b, real sigma, array[] int J, int N, int M, vector sd, vector Z){
        int n_obs = size(Y);
        vector[n_obs] sigma2 = rep_vector(sigma, n_obs);
        return normal_id_glm_marginalized_recover_rng(Y, Xc, mu, b, sigma2, J, N, M, sd, Z);
    }


   /* recovery of the marginalized effect in a normal_id_glm mixed-effect model with varying noise scales where one term with multiple data vectors is marginalized out
   * Args:
   *   Y: a vector of the response variable
   *   Xc: the data matrix in the linear regression part
   *   mu: the mean of the mixed-effect parts that are not marginalized
   *   b: the coefficients in the linear regression part
   *   sigma: a vector of the observation noise scales
   *   J: an array of the grouping indices of the marginalized effect
   *   N: number of groups in the marginalized effect
   *   M: number of data vectors in the marginalized effect, M>1 in this function
   *   sd: a vector of the noise scales for the marginalized effect
   *   L: Cholesky factor of the correlation matrix for the marginalized effect
   *   Z: an array of the data vectors
   * Returns:
   *   a matrix indicating the marginalized effects
   */
    matrix normal_id_glm_multi_marginalized_recover_rng(vector Y, matrix Xc, vector mu, vector b, vector sigma, array[] int J, int N, int M, vector sd, matrix L, array[] vector Z){
        int n_obs = size(Y);
        vector[n_obs] z = Y - mu - Xc * b;
        matrix[M,M] s_u = multiply_lower_tri_self_transpose(diag_pre_multiply(sd, L));
        matrix[M,M] inv_s_u = inverse_spd(s_u);
        vector[n_obs] s_y = sigma .* sigma;
        array[N,M,M] real G;
        for (m1 in 1:M){
            for (m2 in 1:M){
                G[:,m1,m2] = to_array_1d(bincount(Z[m1] .* Z[m2]./ s_y, J));
            }
        }
        array[N,M,M] real F;
        for (n in 1:N){
            F[n] = to_array_2d(to_matrix(G[n]) + inv_s_u);
        }
        array[N,M,M] real F_inv = batch_inverse_spd(F);
        array[N,M] real x;
        for (m in 1:M){
            x[:,m] = to_array_1d(bincount(z .* Z[m] ./ s_y, J));
        }
        matrix[M,N] ans;
        for (n in 1:N){
            ans[:,n] = multi_normal_rng(s_u * (identity_matrix(M) - to_matrix(G[n]) * to_matrix(F_inv[n])) * to_vector(x[n]), to_matrix(F_inv[n]));
        }
        return ans;
    }

   /* recovery of the marginalized effect in a normal_id_glm mixed-effect model with fixed noise scales where one term with multiple data vectors is marginalized out
   * Args:
   *   Y: a vector of the response variable
   *   Xc: the data matrix in the linear regression part
   *   mu: the mean of the mixed-effect parts that are not marginalized
   *   b: the coefficients in the linear regression part
   *   sigma: the fixed observation noise scales
   *   J: an array of the grouping indices of the marginalized effect
   *   N: number of groups in the marginalized effect
   *   M: number of data vectors in the marginalized effect, M>1 in this function
   *   sd: a vector of the noise scales for the marginalized effect
   *   L: Cholesky factor of the correlation matrix for the marginalized effect
   *   Z: an array of the data vectors
   * Returns:
   *   a matrix indicating the marginalized effects
   */
    matrix normal_id_glm_multi_marginalized_recover_rng(vector Y, matrix Xc, vector mu, vector b, real sigma, array[] int J, int N, int M, vector sd, matrix L, array[] vector Z){
        int n_obs = size(Y);
        vector[n_obs] sigma2 = rep_vector(sigma, n_obs);
        return normal_id_glm_multi_marginalized_recover_rng(Y, Xc, mu, b, sigma2, J, N, M, sd, L, Z);
    }
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int<lower=1> Kc;  // number of population-level effects after centering
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  array[N] int<lower=1> J_1;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  vector[N] Z_1_2;
  int<lower=1> NC_1;  // number of group-level correlations
  // data for group-level effects of ID 2
  int<lower=1> N_2;  // number of grouping levels
  int<lower=1> M_2;  // number of coefficients per level
  array[N] int<lower=1> J_2;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_2_1;
  vector[N] Z_2_2;
  int<lower=1> NC_2;  // number of group-level correlations
  // data for group-level effects of ID 3
  int<lower=1> N_3;  // number of grouping levels
  int<lower=1> M_3;  // number of coefficients per level
  array[N] int<lower=1> J_3;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_3_1;
  // data for group-level effects of ID 4
  int<lower=1> N_4;  // number of grouping levels
  int<lower=1> M_4;  // number of coefficients per level
  array[N] int<lower=1> J_4;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_4_1;
  vector[N] Z_4_2;
  int<lower=1> NC_4;  // number of group-level correlations
  // data for group-level effects of ID 5
  int<lower=1> N_5;  // number of grouping levels
  int<lower=1> M_5;  // number of coefficients per level
  array[N] int<lower=1> J_5;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_5_sigma_1;
  // data for group-level effects of ID 6
  int<lower=1> N_6;  // number of grouping levels
  int<lower=1> M_6;  // number of coefficients per level
  array[N] int<lower=1> J_6;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_6_sigma_1;
  // data for group-level effects of ID 7
  int<lower=1> N_7;  // number of grouping levels
  int<lower=1> M_7;  // number of coefficients per level
  array[N] int<lower=1> J_7;  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_7_sigma_1;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] b;  // regression coefficients
  real Intercept;  // temporary intercept for centered predictors
  real Intercept_sigma;  // temporary intercept for centered predictors
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  matrix[M_1, N_1] z_1;  // standardized group-level effects
  cholesky_factor_corr[M_1] L_1;  // cholesky factor of correlation matrix
  vector<lower=0>[M_2] sd_2;  // group-level standard deviations
  matrix[M_2, N_2] z_2;  // standardized group-level effects
  cholesky_factor_corr[M_2] L_2;  // cholesky factor of correlation matrix
  vector<lower=0>[M_3] sd_3;  // group-level standard deviations
  vector<lower=0>[M_4] sd_4;  // group-level standard deviations
  matrix[M_4, N_4] z_4;  // standardized group-level effects
  cholesky_factor_corr[M_4] L_4;  // cholesky factor of correlation matrix
  vector<lower=0>[M_5] sd_5;  // group-level standard deviations
  array[M_5] vector[N_5] z_5;  // standardized group-level effects
  vector<lower=0>[M_6] sd_6;  // group-level standard deviations
  array[M_6] vector[N_6] z_6;  // standardized group-level effects
  vector<lower=0>[M_7] sd_7;  // group-level standard deviations
  array[M_7] vector[N_7] z_7;  // standardized group-level effects
}
transformed parameters {
  matrix[N_1, M_1] r_1;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_1] r_1_1;
  vector[N_1] r_1_2;
  matrix[N_2, M_2] r_2;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_2] r_2_1;
  vector[N_2] r_2_2;
  matrix[N_4, M_4] r_4;  // actual group-level effects
  // using vectors speeds up indexing in loops
  vector[N_4] r_4_1;
  vector[N_4] r_4_2;
  vector[N_5] r_5_sigma_1;  // actual group-level effects
  vector[N_6] r_6_sigma_1;  // actual group-level effects
  vector[N_7] r_7_sigma_1;  // actual group-level effects
  real lprior = 0;  // prior contributions to the log posterior
  // compute actual group-level effects
  r_1 = scale_r_cor(z_1, sd_1, L_1);
  r_1_1 = r_1[, 1];
  r_1_2 = r_1[, 2];
  // compute actual group-level effects
  r_2 = scale_r_cor(z_2, sd_2, L_2);
  r_2_1 = r_2[, 1];
  r_2_2 = r_2[, 2];
  // compute actual group-level effects
  r_4 = scale_r_cor(z_4, sd_4, L_4);
  r_4_1 = r_4[, 1];
  r_4_2 = r_4[, 2];
  r_5_sigma_1 = (sd_5[1] * (z_5[1]));
  r_6_sigma_1 = (sd_6[1] * (z_6[1]));
  r_7_sigma_1 = (sd_7[1] * (z_7[1]));
  lprior += student_t_lpdf(Intercept | 3, 39.5, 46.7);
  lprior += student_t_lpdf(Intercept_sigma | 3, 0, 2.5);
  lprior += cauchy_lpdf(sd_1 | 0,2)
    - 2 * cauchy_lccdf(0 | 0,2);
  lprior += lkj_corr_cholesky_lpdf(L_1 | 1);
  lprior += cauchy_lpdf(sd_2 | 0,2)
    - 2 * cauchy_lccdf(0 | 0,2);
  lprior += lkj_corr_cholesky_lpdf(L_2 | 1);
  lprior += cauchy_lpdf(sd_3 | 0,2)
    - 1 * cauchy_lccdf(0 | 0,2);
  lprior += cauchy_lpdf(sd_4 | 0,2)
    - 2 * cauchy_lccdf(0 | 0,2);
  lprior += lkj_corr_cholesky_lpdf(L_4 | 1);
  lprior += student_t_lpdf(sd_5 | 3, 0, 46.7)
    - 1 * student_t_lccdf(0 | 3, 0, 46.7);
  lprior += student_t_lpdf(sd_6 | 3, 0, 46.7)
    - 1 * student_t_lccdf(0 | 3, 0, 46.7);
  lprior += student_t_lpdf(sd_7 | 3, 0, 46.7)
    - 1 * student_t_lccdf(0 | 3, 0, 46.7);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] sigma = rep_vector(0.0, N);
    mu += Intercept + Xc * b;
    sigma += Intercept_sigma;
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_1_2[J_1[n]] * Z_1_2[n] + r_2_1[J_2[n]] * Z_2_1[n] + r_2_2[J_2[n]] * Z_2_2[n] + r_4_1[J_4[n]] * Z_4_1[n] + r_4_2[J_4[n]] * Z_4_2[n];
    }
    for (n in 1:N) {
      // add more terms to the linear predictor
      sigma[n] += r_5_sigma_1[J_5[n]] * Z_5_sigma_1[n] + r_6_sigma_1[J_6[n]] * Z_6_sigma_1[n] + r_7_sigma_1[J_7[n]] * Z_7_sigma_1[n];
    }
    sigma = exp(sigma);
    target += normal_id_glm_marginalized_lpdf(Y | Xc, mu, b, sigma, J_3, N_3, M_3, sd_3,  Z_3_1);
  }
  // priors including constants
  target += lprior;
  target += std_normal_lpdf(to_vector(z_1));
  target += std_normal_lpdf(to_vector(z_2));
  target += std_normal_lpdf(to_vector(z_4));
  target += std_normal_lpdf(z_5[1]);
  target += std_normal_lpdf(z_6[1]);
  target += std_normal_lpdf(z_7[1]);
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
  // actual population-level intercept
  real b_sigma_Intercept = Intercept_sigma;
  // compute group-level correlations
  corr_matrix[M_1] Cor_1 = multiply_lower_tri_self_transpose(L_1);
  vector<lower=-1,upper=1>[NC_1] cor_1;
  // compute group-level correlations
  corr_matrix[M_2] Cor_2 = multiply_lower_tri_self_transpose(L_2);
  vector<lower=-1,upper=1>[NC_2] cor_2;
  array[M_3] vector[N_3] z_3;  // standardized group-level effects
  // compute group-level correlations
  corr_matrix[M_4] Cor_4 = multiply_lower_tri_self_transpose(L_4);
  vector<lower=-1,upper=1>[NC_4] cor_4;
  // extract upper diagonal of correlation matrix
  for (k in 1:M_1) {
    for (j in 1:(k - 1)) {
      cor_1[choose(k - 1, 2) + j] = Cor_1[j, k];
    }
  }
  // extract upper diagonal of correlation matrix
  for (k in 1:M_2) {
    for (j in 1:(k - 1)) {
      cor_2[choose(k - 1, 2) + j] = Cor_2[j, k];
    }
  }
  // extract upper diagonal of correlation matrix
  for (k in 1:M_4) {
    for (j in 1:(k - 1)) {
      cor_4[choose(k - 1, 2) + j] = Cor_4[j, k];
    }
  }
  // initialize linear predictor term
  vector[N] mu = rep_vector(0.0, N);
  // initialize linear predictor term
  vector[N] sigma = rep_vector(0.0, N);
  mu += Intercept + Xc * b;
  sigma += Intercept_sigma;
  for (n in 1:N) {
    // add more terms to the linear predictor
    mu[n] += r_1_1[J_1[n]] * Z_1_1[n] + r_1_2[J_1[n]] * Z_1_2[n] + r_2_1[J_2[n]] * Z_2_1[n] + r_2_2[J_2[n]] * Z_2_2[n] + r_4_1[J_4[n]] * Z_4_1[n] + r_4_2[J_4[n]] * Z_4_2[n];
  }
  for (n in 1:N) {
    // add more terms to the linear predictor
    sigma[n] += r_5_sigma_1[J_5[n]] * Z_5_sigma_1[n] + r_6_sigma_1[J_6[n]] * Z_6_sigma_1[n] + r_7_sigma_1[J_7[n]] * Z_7_sigma_1[n];
  }
  z_3 = normal_id_glm_marginalized_recover_rng(Y, Xc, mu, b, sigma, J_3, N_3, M_3, sd_3,  Z_3_1);
}
