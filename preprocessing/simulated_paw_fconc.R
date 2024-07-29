
# right now only step functions are incorporated, no intercept
simulate_fconc = function(I = 50, D = 100, 
                          sigma = 0.1, 
                          seed = 546474, 
													t_max = 1, 
													rand_int = FALSE, 
													lambda0 = 40){
	
	set.seed(seed)
	Kt = 10
	
	# define grid and basis
	# intercept only included if x = "periodic"
	grid = seq(0, t_max, length.out = D)
	
	
	# upper triangle will be zeroed out
	beta_coefs_x = .1 * c(2, 2, 30, 28, 25, 18, 15, 10, 7, 7)
	beta_coefs_0 = c(0,2, 3, 7, 6, 1, 5, 0, 1, 0)

	
	beta_basis = bs(grid, df = Kt, intercept = TRUE)
	beta_0 = beta_basis %*% beta_coefs_0
	beta_x = beta_basis %*% beta_coefs_x
	
	xmat = matrix(NA, nrow = I, ncol = D)
	jump_points = runif(I, quantile(grid, .1), quantile(grid, .9))
	#sine_shift = runif(I, -0.3, 0.2)
	sine_shift = runif(I, -0.8, 1.4)
	sine_scale = rnorm(I, 1, 0.5)
	sqrt_scale = runif(I, 2, 5)

	# define other elements
	epsilon = matrix(rnorm(I*D, 0, sqrt(sigma)), nrow = I, ncol = D)
	Y = matrix(NA, I, D)
	colnames(Y) = paste0("time_", 1:D)
	
	lambda = rep(lambda0, Kt * I)
	re_coef = matrix(rnorm(I * Kt, 0, sqrt(lambda)), nrow = I, ncol = Kt)
	re_mat = re_coef %*% t(beta_basis)
	re_mat_data = matrix(NA, I, D)	
	xstar_beta = matrix(NA, I, D)
	
	
	for(i in 1:I){
		xmat[i,] = sine_scale[i] * sin(pi*seq(0, 1, length.out = D) + 1 + sine_shift[i])
		xstar_beta[i, ] =  beta_x * xmat[i,]
		
		Y[i,] = beta_0 +
			xstar_beta[i, ] +
			epsilon[i,]	
		
		if(rand_int){
			re_mat_data[i, ] = re_mat[i,]
			
			Y[i,] = beta_0 +
				xstar_beta[i, ] +
				re_mat_data[i, ] +
				epsilon[i,]	
		}
		
	} # end for loop
	
	
	simulated_data = data.frame(
		trial = rep(1:I, each = D),
		time = rep(grid, I),
		value = as.vector(t(Y)),
		x = as.vector(t(xmat)),
		re = as.vector(t(re_mat_data)),
		xBeta = as.vector(t(xstar_beta))
	)
	
	coef_df = data.frame(
		time = grid,
		beta0 = beta_0,
		beta_1 = beta_x,
		trial = 0
	)
	
	list(data = simulated_data, 
	     epsilon = epsilon,
			 coefficient_fns = coef_df, 
			 surface = NULL, 
			 re_mat_data = re_mat_data)
	
} # end function

# This function integrates the coefficient surface
integrate_surface = function(basis_time, beta_f, covar = NULL){
  D = length(basis_time)
  if(is.null(covar)){
    covar = rep(1, length.out = D)
  }
  
  beta_covar_f = function(s, t = 1, basis_time, beta_f, covar){
     beta_f(s,t) *  (s < t) *  covar[which(basis_time == s)] 
  }
  
  Gstar = sapply(basis_time, FUN = function(t){cumtrapz(basis_time, beta_covar_f(basis_time, t, basis_time, beta_f,covar))[D,]})
  return(t(Gstar))
}

# returns a matrix which is a surface where the columns are S and the rows are T
# equals zero when s > t
make_sim_surface = function(beta_f, D, t_max){
	t = s = seq(0, t_max, length.out = D)
	surface = matrix(NA, D, D)
	colnames(surface) = paste0("s", s)
	rownames(surface) = paste0("t", t)
	
	for(time_index in 1:D){
		surface[,time_index] = beta_f(s[time_index], t) * (s[time_index] <= t)
	}
	
	surface
}




x_f = function(grid, jump_point){
	ifelse(grid < jump_point, 0, 1)
}


