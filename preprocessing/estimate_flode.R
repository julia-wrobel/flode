
# forcing functions are a design matrix, if you want to include an intercept it needs to include one
estimate_flode = function(data, Kt = 10, 
													alpha0 = 6,  
													forcing_functions = c(""), 
													max_iter = 100,
													tol =  0.00001,
													sigma_b = 100, 
													sigma_d = 100,
													random_effects = TRUE,
													initial_position = TRUE,
													y0_vec = NULL,
													a = 0.001,
													alpha_upper = 14){
	# process data to extract forcing functions
	x_data = data %>% select(trial, all_of(forcing_functions))  %>% nest(data = c(-trial)) %>% pull(data)
	
	basis_time = filter(data, trial == first(trial))$time # assumes grid is even across subjects
	spline_basis = bs(basis_time, df = Kt, intercept = TRUE)
	
	N = length(unique(data$trial))
	D = length(basis_time)
	P = length(forcing_functions)
	trial_ids = unique(data$trial) 
	
	##################################################
	# initialize variance parameters 
	sigma_cur = 1
	sigma_d_cur = sigma_d
	sigma_b_cur = sigma_b
	
	##################################################
	# define penalty matrix
	diff2 = matrix(rep(c(1,-2,1, rep(0, D-2)), D-2)[1:((D-2)*D)], D-2, D, byrow = TRUE)
	P2 = t(spline_basis) %*% t(diff2) %*% diff2 %*% spline_basis # not full rank
	Pen = a * diag(Kt) + (1 - a) * P2
	
	# make yiStar term
	if(is.null(y0_vec)){
		y0_vec = (data %>% group_by(trial) %>% slice(1) %>% ungroup())$value # data should be sorted by time for each subject
	}else{
		y0_vec = y0_vec
	}
	y0_cur = rep(y0_vec, each = D)
	data = data %>% mutate(Yo = y0_cur, 
												 y0_star = Yo * exp(- alpha0 * time))
	
	##################################################
	# initialize alpha and beta
	alpha_cur = alpha0
	thetaStar = integrate_alpha(basis_time, spline_basis, alpha_cur, covar = NULL)
	xstar_ls = lapply(x_data, make_Xstar, alpha = alpha_cur, spline_basis = spline_basis, basis_time = basis_time,
										thetaStar = thetaStar)
	xstar = do.call(rbind, xstar_ls)
	beta_coefs = solve(crossprod(xstar)) %*% (t(xstar) %*% (data$value - data$y0_star))
	beta_coefs = matrix(beta_coefs, nrow = Kt, ncol = P)
	beta_cur = spline_basis %*% beta_coefs
	
	##################################################
	# initialize coefficients for random intercepts
	#delta_coefs = matrix(rnorm(Kt *  N, 0, sqrt(sigma_d_cur)), nrow = Kt, ncol = N)
	delta_coefs = matrix(0, nrow = Kt, ncol = N)
	deltaP_squared = matrix(0, nrow = 1, ncol = N)
	Dstar = thetaStar
	C = matrix(0, nrow = Kt, ncol = Kt)
	data = data %>% mutate(delta =  as.vector(spline_basis %*% delta_coefs),
												 deltaStar = as.vector(Dstar %*% delta_coefs))
	
	
	##################################################
	# iteratively calculate alpha and beta
	iter = 1
	error_cur = 1000
	sse_vec = alpha_loss_vec = logLik_vec = alpha_vec = sigma_b_vec = sigma_d_vec = sigma_vec = rep(NA, max_iter)
	beta0_mat = beta1_mat = matrix(NA, nrow = max_iter, ncol = Kt)
	
	## save initial parameter values
	alpha_vec[1] = alpha_cur
	sigma_b_vec[1] = sigma_b_cur
	sigma_d_vec[1] = sigma_d_cur
	sigma_vec[1] = sigma_cur
	# not updating C yet		
	beta0_mat[1,] = beta_coefs[, 1]
	beta1_mat[1,] = beta_coefs[, 2]
	sse_direction = NA
	
	while(iter < max_iter && abs(error_cur) > tol){
	#while(iter < max_iter ){
		##################################################
		# update alpha
		alpha_nls = optim(alpha_cur, 
											alpha_loss_re, 
											data = data, 
											x_data = x_data,
											spline_basis = spline_basis, 
											basis_time = basis_time,
											beta_coefs = beta_coefs, 
											delta_coefs = delta_coefs, 
											random_effects = random_effects,
											C = C, # part of random effects but shows up in error term
											method = "Brent", 
											lower = max(0, alpha_cur-1), upper = alpha_cur+1)
		alpha_loss_vec[iter] = alpha_nls$value
		alpha_diff = abs(alpha_cur - alpha_nls$par)
		alpha_cur = alpha_nls$par 
			
		##################################################
		# update sigma
		RSS = alpha_nls$value
		sigma_cur = as.numeric(RSS)/(N * D)
		
		##################################################
		# update bases containing alpha 
		# these are all the numeric integrals. Ultimately we may want to reparameterize and calculate these in a more efficient way
		data = data %>% mutate(y0_star = Yo * exp(- alpha_cur * time))
		thetaStar = integrate_alpha(basis_time, spline_basis, alpha_cur, covar = NULL)
		xstar_ls = lapply(x_data, make_Xstar, alpha = alpha_cur, spline_basis = spline_basis, 
											basis_time = basis_time,
											thetaStar = thetaStar)
		xstar = do.call(rbind, xstar_ls)
		
		if(random_effects & iter > 2){	
			Dstar = thetaStar
			C = solve(1/sigma_d_cur * Pen + crossprod(Dstar)/sigma_cur) 
			trace_PC = sum(diag(Pen %*% C))
			C_Dstar = C %*% t(Dstar)
		}else{
			C = matrix(0, nrow = Kt, ncol = Kt)
		}
		
		
		##################################################
		# update beta
		beta_coefs = solve(crossprod(xstar) + sigma_cur / (2 * sigma_b_cur) * Matrix::bdiag(rep(list(Pen + t((Pen))), P))) %*% (t(xstar) %*% (data$value - data$y0_star - as.vector(Dstar %*% delta_coefs)))
		beta_coefs  = matrix(beta_coefs, nrow = Kt, ncol = P)
		beta_cur = spline_basis %*% beta_coefs
		
		##################################################
		# update sigma_b
		sigma_b_cur = sum(diag(crossprod(beta_coefs, Pen %*% beta_coefs)))/(Kt*P)

		##################################################
		# update initial position
		
		if(initial_position){
			e_alpha_t = exp(-alpha_cur * basis_time)
			e_2alpha_t = sum(exp(- 2 * alpha_cur * basis_time))
			
			for(i in 1:length(trial_ids)){
				trial_df = filter(data, trial == trial_ids[i])
				
				y0_vec[i] = sum(e_alpha_t * (trial_df$value - trial_df$deltaStar - xstar_ls[[i]] %*% as.vector(beta_coefs))) / e_2alpha_t
			}
			y0_cur = rep(y0_vec, each = D)
			
			
			data = data %>% mutate(Yo = y0_cur,
														 y0_star = Yo * exp(- alpha_cur * time))
		}
		
		####################################################################################################
		####################################################################################################
		####################################################################################################
		# E-step (see E-step of paper). Update delta and expected values containing delta
		####################################################################################################
		####################################################################################################
		####################################################################################################
		
		# Define global updates outside of loop to avoid recomputing each time
		if(random_effects & iter > 2){
			# loop over trials
			for(i in 1:length(trial_ids)){
				trial_df = filter(data, trial == trial_ids[i])
				
				delta_coefs[, i] =   (C_Dstar %*% (trial_df$value - trial_df$y0_star -
																					 	xstar_ls[[i]] %*% as.vector(beta_coefs)))/sigma_cur
				
				deltaP_squared[, i] = trace_PC + crossprod(delta_coefs[,i], Pen %*% delta_coefs[,i])
			}
			
			
			## update terms that contain delta coefs
			sigma_d_cur = sum(deltaP_squared)/(N * Kt)
		}
		
			
		
		# probably not necessary to calculate yhat here
		# remove if decided not necessary
		data = data %>% mutate(delta = as.vector(spline_basis %*% delta_coefs),
													 deltaStar = as.vector(Dstar %*% delta_coefs),
													 yhat = y0_star + as.numeric(xstar %*% as.vector(beta_coefs)) + deltaStar)
		
		sse_vec[iter] = sum((data$value - data$yhat)^2)
		
		# calculate convergence parameters
		##################################################
		logLikelihood = logLik(data = data, N = N, D = D, Kt = Kt,
													 Dstar = Dstar, delta_coefs = delta_coefs,
													 beta_coefs = beta_coefs, 
													 xstar_ls = xstar_ls, 
													 sigma = sigma_cur, 
													 sigma_d = sigma_d_cur, Pen = Pen,
													 random_effects = random_effects)

		
		logLik_vec[iter] = (logLikelihood$ylik)/(N * D)
		error_cur = ifelse(iter < 2, 1000, diff(logLik_vec)[iter - 1])
		
		
		### remove this and replace with thing from bayes_fosr
		message(paste0("current iteration: ", iter))
		message(paste0("current alpha: ", alpha_cur))
		#message(paste0("current likelihood error: ", error_cur))
		#message(paste0("current log Likelihood: ", logLik_vec[iter]))
		
		
		############################################################
		## save all parameter values
		alpha_vec[iter + 1] = alpha_cur
		sigma_b_vec[iter + 1] = sigma_b_cur
		sigma_d_vec[iter + 1] = sigma_d_cur
		sigma_vec[iter + 1] = sigma_cur
	
		beta0_mat[iter + 1,] = beta_coefs[, 1]
		beta1_mat[iter + 1,] = beta_coefs[, 2]
		
		############################################################
		#
		# 
		iter = iter + 1

		if(error_cur < 0 & iter > 55){
			break
		}
		
		if(iter > 55){
			sse_direction = sign(sse_vec[iter-1] - sse_vec[iter-2])
			if(sse_direction == 1){
				break
			}
		}

	} # end while loop
	##################################################
	
	data = data %>% mutate(yhat = y0_star + as.numeric(xstar %*% as.vector(beta_coefs)) + deltaStar)

	
	############################################################
	# get surface for each forcing function
	surfaces = as.list(rep(NA, P))
	for(p in 1:P){
		beta_p = beta_cur[, p]
		surfaces[[p]] = make_surface(basis_time, alpha_cur, beta = beta_p) 
	}
	
	colnames(beta_cur) = paste0("beta", 0:(P-1))
	beta_cur = as_tibble(beta_cur) %>% mutate(time = basis_time)
	
	
	
	# variance of beta coefficients
	G = Matrix::bdiag(rep(list(sigma_b_cur * solve(Pen)), P)) # G is actually symmetric, at least it should be
	A = solve(crossprod(xstar) + sigma_cur * solve(G)) %*% t(xstar)
	var_y_given_d = xstar %*% G %*% t(xstar) + sigma_cur * diag(N * D)
	var_b_given_d = A %*% var_y_given_d %*% t(A) + G - 2*A %*% xstar %*% G
	
	#var_y = sigma_cur*diag(D) +sigma_d_cur * Dstar %*% solve(Pen) %*% t(Dstar)
	#var_y = Matrix::bdiag(rep(list(var_y), N))
	#var_b = A %*% var_y %*% t(A)
	var_components = list(
	  var_b_given_d = var_b_given_d,
	  xstar = xstar,
	  Pen = Pen,
	  P = P,
	  N = N,
	  D = D,
	  #spline_basis = spline_basis,
	  sigma = sigma_cur,
	  alpha = as.numeric(alpha_cur),
	  beta = beta_cur,
	  sigma_d = sigma_d_cur,
	  sigma_b = sigma_b_cur
	)
	
	
	list(data = data, 
			 beta = beta_cur, 
			 beta_coefs = beta_coefs,
			 xB = as.numeric(xstar %*% as.vector(beta_coefs)),
			 alpha = as.numeric(alpha_cur),
			 sigma = sigma_cur,
			 sigma_d = sigma_d_cur,
			 sigma_b = sigma_b_cur,
			 surface = surfaces, 
			 maxiter = iter,
			 y0 = y0_vec,
			 sse = sse_vec[iter - 1],
			 spline_basis = spline_basis,
			 logLikelihood = logLikelihood$ylik,
			 var_components = var_components
			 )
}
