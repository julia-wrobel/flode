library(mvtnorm)
# right now only step functions are incorporated, no intercept
simulate_fhist = function(I = 50, D = 100, 
                          sigma = 0.1, 
                          seed = 546474, 
													t_max = 1, 
													rand_int = FALSE, 
													lambda0 = 40){
	
	set.seed(seed)
	Kt = 10
	
	# define grid and basis
	# intercept only included if x = "periodic"
	#grid = seq(0, t_max, length.out = D)
	grid = 1:D
	
	# upper triangle will be zeroed out
	#beta_coefs_x = t(1:10 *t(matrix(rep(c(2, 2, 30, 28, 25, 18, 15, 10, 7, 7), 10), ncol = 10, nrow = 10)))
	beta_coefs_0 = c(0,2, 3, 7, 6, 1, 5, 0, 1, 0)

	
	beta_basis = bs(grid, df = Kt, intercept = TRUE)
	beta_0 = beta_basis %*% beta_coefs_0
	#beta_x = beta_basis %*% beta_coefs_x
	
	# get P and forcing functions
	P = 2
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
	
	surface_df =
	  expand_grid(
	    s = seq(0, 1, length = D),
	    t = seq(0, 1, length = D)
	  ) |> 
	  mutate(
	    coef = map2_dbl(s, t, \(x,y) dmvnorm(x = c(x, y), mean = c(.25, .75), sigma = matrix(c(.0075, 0, 0, .0075), 2, 2))),
	    coef = ifelse(s < t, coef, NA)
	  )
	
	lambda = rep(lambda0, Kt * I)
	re_coef = matrix(rnorm(I * Kt, 0, sqrt(lambda)), nrow = I, ncol = Kt)
	re_mat = re_coef %*% t(beta_basis)
	re_mat_data = matrix(NA, I, D)	
	xstar_beta = matrix(NA, I, D)
	
	
	for(i in 1:I){
		xmat[i,] = sine_scale[i] * sin(pi*seq(0, 1, length.out = D) + 1 + sine_shift[i])
		xstar_beta[i, ] = integrate_surface(grid, surface_df, covar = xmat[i,]) 
		
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
		trial = 0
	)
	
	list(data = simulated_data, 
	     epsilon = epsilon,
			 coefficient_fns = coef_df, 
			 surface = surface_df, 
			 re_mat_data = re_mat_data)
	
} # end function



integrate_surface = function(basis_time, surface_df, covar){
  D = length(basis_time)
  
  xb_mat = matrix(NA, nrow = D, ncol = D)
  s_grid = unique(surface_df$s)
  for(ss in 1:D){
    xb_mat[,ss] = filter(surface_df, s == s_grid[ss])$coef * covar[ss]
  }
  dif = diff(basis_time)[1]
  return(colSums(xb_mat * dif, na.rm = TRUE))
}




