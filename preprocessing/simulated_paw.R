
# right now only step functions are incorporated, no intercept
simulate_flode = function(I = 50, D = 100, sigma_y0 = 5, sigma = 0.1, 
                          seed = 546474, 
                          alpha = 2,
													t_max = 1, 
													rand_int = FALSE, 
													lambda0 = 40, 
													autoregressive = FALSE,
													constant_beta = FALSE){
	
	set.seed(seed)
	Kt = 10
	
	# define grid 
	grid = seq(0, t_max, length.out = D)
	
	# define basis coefficients
	if(constant_beta){
	  beta_coefs_x = rep(10,10)
	  beta_coefs_0 = rep(0,10)
	  #beta_coefs_0 = c(-5, 5, seq(10, 45, length.out = 5),  40, 40, 40)
	}else{
	  beta_coefs_x = 1.5 * c(2, 2, 30, 28, 25, 18, 15, 10, 7, 7)
	  beta_coefs_0 = c(-5, 5, seq(10, 45, length.out = 5),  40, 40, 40)
	  #beta_coefs_0 = rep(0,10)
	}
	
	# define beta and basis
	beta_basis = bs(grid, df = Kt, intercept = TRUE)
	beta_0 = beta_basis %*% beta_coefs_0
	beta_x = beta_basis %*% beta_coefs_x
	
	# get P (number of forcing functions including intercept) and forcing functions
	P = 2
	xmat = matrix(NA, nrow = I, ncol = D)
	jump_points = runif(I, quantile(grid, .1), quantile(grid, .9))
	sine_shift = runif(I, -0.8, 1.4)
	sine_scale = rnorm(I, 1, 0.5)
	sqrt_scale = runif(I, 2, 5)

	# define other elements
	# Y0 is 
	Y0 = rnorm(I, 0, sqrt(sigma_y0))
	epsilon = matrix(rnorm(I*D, 0, sqrt(sigma)), nrow = I, ncol = D)
	if(autoregressive){
	  epsilon = matrix(NA, nrow = I, ncol = D)
	  for(i in 1:I){
	    epsilon[i,] = arima.sim(list(order=c(1,0,0), ar=.5, sd = sqrt(sigma)), n = D)
	  }
	}
	Y = matrix(NA, I, D)
	colnames(Y) = paste0("time_", 1:D)
	
	### define surface for intercept and beta1
	intercept_surface = make_sim_surface(alpha = alpha, beta = beta_0, D = D, t_max = t_max)
	surface = make_sim_surface(alpha = alpha, beta = beta_x, D = D, t_max = t_max)
	###
	
	# lambda controls variance of random effects
	lambda = rep(lambda0, Kt * I)
	
	# define random effects coefficients
	re_coef = matrix(rnorm(I * Kt, 0, sqrt(lambda)), nrow = I, ncol = Kt)
	re_mat = re_coef %*% t(beta_basis)
	
	# set up matrices to store random effects and xB
	re_mat_data = matrix(NA, I, D)	
	xstar_beta = matrix(NA, I, D)
	
	thetaStar = integrate_alpha(grid, beta_basis, alpha)
	xstar_beta0 = thetaStar %*% beta_coefs_0
	
	for(i in 1:I){
	  # define X
		xmat[i,] = sine_scale[i] * sin(pi*seq(0, 1, length.out = D) + 1 + sine_shift[i])

		# integrate over X * beta
		xstar_beta[i, ] = integrate_alpha(grid, beta_basis, alpha, covar = xmat[i,]) %*% beta_coefs_x
		
		# calculate Y values
		Y[i,] = Y0[i] * exp(-alpha * grid) + 
			xstar_beta0 +
			xstar_beta[i, ] +
			epsilon[i,]	
		
		if(rand_int){
			re_mat_data[i, ] = thetaStar %*% re_coef[i,]
			
			Y[i,] = Y0[i] * exp(-alpha * grid) + 
				xstar_beta0 +
				xstar_beta[i, ] +
				re_mat_data[i, ] +
				epsilon[i,]	
		}
		
	} # end for loop
	
	# save  data results
	simulated_data = data.frame(
		trial = rep(1:I, each = D),
		y0 = rep(Y0, each = D),
		time = rep(grid, I),
		value = as.vector(t(Y)),
		x = as.vector(t(xmat)),
		re = as.vector(t(re_mat)),
		reStar = as.vector(t(re_mat_data)),
		xStarbeta = as.vector(t(xstar_beta))
	)
	
	# save true coefficient functions
	coef_df = data.frame(
		time = grid,
		beta0 = beta_0,
		beta1 = beta_x,
		trial = 0
	)
	
	list(data = simulated_data, alpha = alpha, epsilon = epsilon,
			 y0 = Y0, coefficient_fns = coef_df, surface = surface, 
			 re_mat_data = re_mat_data)
	
} # end function


# returns a matrix which is a surface where the columns are S and the rows are T
# equals zero when s > t
make_sim_surface = function(alpha, beta, D, t_max){
	t = s = seq(0, t_max, length.out = D)
	surface = matrix(NA, D, D)
	colnames(surface) = paste0("s", s)
	rownames(surface) = paste0("t", t)
	
	for(time_index in 1:D){
		surface[,time_index] = exp(-alpha * (t-s[time_index]) ) * beta[time_index] * (s[time_index] <= t)
	}
	
	surface
}

x_f = function(grid, jump_point){
	ifelse(grid < jump_point, 0, 1)
}


