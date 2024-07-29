# used to generate beta surface that is used for comparison with functional regression
make_surface = function(basis_time, alpha, beta_s = 1){
	###
	exp_f = function(s, t, alpha){
		exp(-alpha * (t - s)) * (s <= t)
	}
	
	D = length(basis_time)
	s = t = basis_time
	
	exp_surface = t(outer(s, t, exp_f, alpha = alpha) * beta_s)
	colnames(exp_surface) = paste0("s", 1:D); rownames(exp_surface) = paste0("t", 1:D)
	
	attr(exp_surface, "t") = t; attr(exp_surface, "s") = s
	return(exp_surface)
}

## function to get joint log likelihood of model for monitoring convergence
logLik = function(data, N, D, Kt, Dstar, delta_coefs, beta_coefs, xstar_ls, sigma, sigma_d, Pen,
									random_effects){
	
	if(random_effects){
		lik_y_sigma = sigma * diag(D)  + Dstar %*% (sigma_d * solve(Pen)) %*% t(Dstar)
	}else{
		lik_y_sigma = sigma * diag(D)  
	}
	
	lik_y = rep(NA, length.out = N)
	# loop over trials
	trial_ids = unique(data$trial)
	for(i in 1:length(trial_ids)){
		trial_df = filter(data, trial == trial_ids[i])
		mu = (trial_df$y0_star + as.numeric(xstar_ls[[i]] %*% as.vector(beta_coefs)))
		
		lik_y[i] = dmvnorm(trial_df$value, mean = mu, sigma = lik_y_sigma,
											 checkSymmetry = FALSE,
											 log = TRUE) 
	}
	
	list(ylik = sum(lik_y))

}


# this function puts together Xstar for ith trial
# thetaStar is stored value of integrated intercept term at a particular value of alpha
make_Xstar = function(xdata, alpha, spline_basis, basis_time, thetaStar){
	forcing_functions = as.matrix(xdata)
	D = nrow(forcing_functions)
	P = ncol(forcing_functions) # includes intercept term
	Kt = ncol(spline_basis)
	
	xstar_i = matrix(NA, nrow = D, ncol = P * Kt)
	for(p in 1:P){
		p_index = (p-1) * Kt + 1
		
		if(p > 1){
			xstar_ip = integrate_alpha(basis_time, spline_basis, alpha, covar = forcing_functions[,p])
		}else{
			xstar_ip = thetaStar
		}
		xstar_i[, p_index:(p * Kt)] = xstar_ip
	} # end for loop
	
	return(xstar_i)
}

# This function integrates the exponential surface containing alpha
# used to create Dstar and Xstar terms
# returns Gstar, a D x Kt matrix that will be multiplied by beta_coefs or delta_coefs in later steps
# covar is either a length D vector of ones or covariate values for a particular covariate and subjects
integrate_alpha = function(basis_time, spline_basis, alpha, covar = NULL){
	D = length(basis_time)
	if(is.null(covar)){
		covar = rep(1, length.out = D)
	}
	
	exp_f = function(s, t = 1, basis_time, alpha, covar){
		exp(-alpha * (t - s))  *  (s < t) *  covar[which(basis_time == s)] * spline_basis[which(basis_time == s),]
	}
	
	Gstar = sapply(basis_time, FUN = function(t){cumtrapz(basis_time, exp_f(basis_time, t, basis_time, alpha, covar))[D,]})
	return(t(Gstar))
}


# alpha loss for no random effects
alpha_loss = function(data, x_data, beta_coefs, alpha, spline_basis){
	
	data = data %>% mutate(y0_star = Yo * exp(- alpha * time))
	
	basis_time = filter(data, trial == 1)$time
	thetaStar = integrate_alpha(basis_time, spline_basis, alpha, covar = NULL)
	xstar_ls = lapply(x_data, make_Xstar, alpha = alpha, spline_basis = spline_basis, basis_time = basis_time, thetaStar = thetaStar)
	xstar = do.call(rbind, xstar_ls)
	
	#Yo * exp(- alpha0 * time)
	#unique(subject_df$Yo) * exp(-alpha * subject_df$time)
	
	sum((data$value - data$y0_star - as.numeric(xstar %*% as.vector(beta_coefs)))^2)
}


alpha_loss_re = function(data, x_data, spline_basis, basis_time, beta_coefs, alpha, 
												 random_effects, delta_coefs, C){
	
	if(!random_effects){
		alpha_loss(data, x_data, beta_coefs, alpha, spline_basis)
	}else{
		Kt = dim(spline_basis)[2]
		Dstar = integrate_alpha(basis_time, spline_basis, alpha, covar = NULL) # same as thetaStar
		deltaStar = Dstar %*% delta_coefs
		Dstar_squared = crossprod(Dstar)
		trace_Dstar_C = sum(diag(Dstar_squared %*% C))
		trial_ids = unique(data$trial)
		xstar_ls = lapply(x_data, make_Xstar, alpha = alpha, 
											spline_basis = spline_basis, basis_time = basis_time, 
											thetaStar = Dstar)
		#xstar = do.call(rbind, xstar_ls)
		
		#sum((data$value - data$y0_star - as.numeric(xstar %*% as.vector(beta_coefs)) - deltaStar)^2)
		
		rss = rep(NA, length(trial_ids))
		for(i in 1:length(trial_ids)){
			subject_df = filter(data, trial == trial_ids[i])
			ystar_i = unique(subject_df$Yo) * exp(-alpha * subject_df$time)
			deltastar_i = deltaStar[, i]
			# maybe below shouldn't be subject specific or maybe tehre is a plus when there should be a minus
			deltaStar_squared_i = trace_Dstar_C + crossprod(delta_coefs[,i], Dstar_squared %*% delta_coefs[,i])
			Y_centered = subject_df$value - ystar_i - as.numeric(xstar_ls[[i]]  %*% as.vector(beta_coefs))


			rss[i] = crossprod(Y_centered) -
				2 * crossprod(Y_centered, deltastar_i) +
				deltaStar_squared_i

		}
		return(sum(rss))
	}
	
	
}



# used to generate beta surface that is used for comparison with functional regression
get_ests = function(f){
  load(f)
  
  df = as_tibble(results$params) %>%
    mutate(
      time_flode = results$times[[1]]$toc -
        results$times[[1]]$tic,
      time_fhist = results$times[[2]]$toc -
        results$times[[2]]$tic,
      time_fconc = results$times[[3]]$toc -
        results$times[[3]]$tic
    )
  
  if(!grepl("genFhist", f)){
    simulated_data = simulate_flode(I = df$N, 
                                    D = df$D, 
                                    sigma_y0 = 5, 
                                    sigma = 0.1, 
                                    alpha = df$alpha_true, 
                                    rand_int = TRUE, 
                                    lambda0 = 20,
                                    seed = df$seed)
    df$genMod = "flode"
    df$ISEb0 = mean(results$flode_results$beta$beta0 - simulated_data$coefficient_fns$beta0)
    df$ISEb1 = mean(results$flode_results$beta$beta1 - simulated_data$coefficient_fns$beta1)
    true_surface = simulated_data$surface
    
    # calculate surface errors
    fhist_surface = coef(results$fhist_results, n1 = df$D, n2 = df$D)$smterms$"ff(x,s)"$coef %>%
      select(s = x.smat, t = x.tmat, value = value) %>% as_tibble() 
    
    
  }else{
    source(here::here("preprocessing",
                      "pffr_sim.R"))
    set.seed(df$seed)
    simulated_data = pffr_sim(scenario = c("int", "ff"),
                              df$N, nxgrid = df$D, SNR = 10,
                              sigma_d = 1,
                              rand_int = TRUE)
    
    df$genMod = "fhist"
    true_surface = t(simulated_data$surface)
    
    # calculate surface errors
    fhist_surface = coef(results$fhist_results, n1 = df$D, n2 = df$D)$smterms$"ff(X1,s)"$coef %>%
      select(s = X1.smat, t = X1.tmat, value = value) %>% as_tibble() 
  }

  
  ## interpolate surface to get on same grid as flode
  interp_obj = list(x = unique(fhist_surface$s), y = unique(fhist_surface$t),
                    z = matrix(fhist_surface$value, df$D, df$D, byrow = TRUE))
  tt = ss = seq(0, 1, length.out = df$D)
  loc_mat = list(x = ss, y = tt)
  fhist_surface = fields::interp.surface.grid(interp_obj, loc_mat)$z
  fhist_surface[df$D,] = fhist_surface[df$D-1,] # carry last value forward so not NA
  fhist_surface = fhist_surface * lower.tri(fhist_surface, diag = TRUE)# zero out upper triangle of the matrix
  
  fhist_intercept = coef(results$fhist_results, n1 = df$D, n2 = df$D)$smterms$"Intercept(t)"$coef %>%
    select(t = t.vec, value = value) %>% as_tibble() 
  
  
  df = df %>%
    mutate(alpha_hat = results$flode_results$alpha,
           surfaceErr_flode =  sum((true_surface -                           
                                      as.vector(results$flode_results$surface[[2]]))^2,
                                   na.rm = TRUE),
           surfaceErr_fhist = sum((true_surface - fhist_surface)^2,
                                  na.rm = TRUE),
           sse_flode = results$flode_results$sse,
           sse_fhist = sum((results$fhist_results$residuals)^2),
           sse_fconc =  sum((results$fconc_results$data$value -
                               results$fconc_results$Yhat)^2))
  
  df
  # end function  
}  



###
# process files for simulation results
get_ests_fhist_only = function(f){
  load(f)
  
  df = as_tibble(results$params) %>%
    mutate(
      time_fhist = results$times[[1]]$toc -
        results$times[[1]]$tic,
    )
  
  simulated_data = simulate_flode(I = df$N, 
                                  D = df$D, 
                                  sigma_y0 = 5, 
                                  sigma = 0.1, 
                                  alpha = df$alpha_true, 
                                  rand_int = TRUE, 
                                  lambda0 = 20,
                                  seed = df$seed)
  df$genMod = "flode"
  true_surface = simulated_data$surface
  
  # calculate surface errors
  fhist_surface = coef(results$fhist_results, n1 = df$D, n2 = df$D)$smterms$"ff(x,s)"$coef %>%
    select(s = x.smat, t = x.tmat, value = value) %>% as_tibble() 
  
  ## interpolate surface to get on same grid as flode
  interp_obj = list(x = unique(fhist_surface$s), y = unique(fhist_surface$t),
                    z = matrix(fhist_surface$value, df$D, df$D, byrow = TRUE))
  tt = ss = seq(0, 1, length.out = df$D)
  loc_mat = list(x = ss, y = tt)
  fhist_surface = fields::interp.surface.grid(interp_obj, loc_mat)$z
  fhist_surface[df$D,] = fhist_surface[df$D-1,] # carry last value forward so not NA
  fhist_surface = fhist_surface * lower.tri(fhist_surface, diag = TRUE)# zero out upper triangle of the matrix

  # colnames(fhist_surface) = paste0("s", ss)
  # fhist_surface = fhist_surface %>%
  #   as_tibble() %>%
  #   mutate(t = tt) %>% select(t, everything()) %>%
  #   gather(s, value, starts_with("s")) %>%
  #   arrange(t) %>%
  #   mutate(s = rep(ss, df$D),
  #          alpha = df$alpha_true,
  #          method = "fhist") %>%
  #   mutate(value = ifelse(s <= t, value, NA))

  
  fhist_intercept = coef(results$fhist_results, n1 = df$D, n2 = df$D)$smterms$"Intercept(t)"$coef %>%
    select(t = t.vec, value = value) %>% as_tibble() 
  
  
  df = df %>%
    mutate(surfaceErr_fhist = sum((true_surface - fhist_surface)^2,
                                  na.rm = TRUE),
           sse_fhist = sum((results$fhist_results$residuals)^2))
  
  df
  # end function 
}



get_surfaces = function(f, f_fhist){
  load(f_fhist)
  
  # get true surface
  df = as_tibble(results$params)
  tt = ss = seq(0, 1, length.out = df$D)
  true_surface = results$surface %>% as_tibble() %>%
    mutate(t = tt) %>% select(t, everything()) %>%
    gather(s, value, starts_with("s")) %>%
    arrange(t) %>%
    mutate(s = rep(ss, df$D),
           value = ifelse(s <= t, value, NA)) %>%  
    mutate(alpha = df$alpha_true,
           method = "truth")
  
  ## calculate surface for fhist
  fhist_surface = coef(results$fhist_results, n1 = df$D, n2 = df$D)$smterms$"ff(x,s)"$coef %>%
    select(s = x.smat, t = x.tmat, value = value) %>% 
    as_tibble()  %>%
    #mutate(s = round(s, 4)) #%>%
  mutate(#value = ifelse(s >= t, NA, value),
         alpha = df$alpha_true,
         method = "fhist") %>%
  select(s,t, value, alpha, method)
  
  ## interpolate surface to get on same grid as flode
  interp_obj = list(x = unique(fhist_surface$s), y = unique(fhist_surface$t),
                    z = matrix(fhist_surface$value, df$D, df$D, byrow = TRUE))
  
  loc_mat = list(x = ss, y = tt)
  fhist_surface = fields::interp.surface.grid(interp_obj, loc_mat)$z
  fhist_surface[df$D,] = fhist_surface[df$D-1,] # carry last value forward so not NA
  fhist_surface = fhist_surface * lower.tri(fhist_surface, diag = TRUE)# zero out upper triangle of the matrix
  
  colnames(fhist_surface) = paste0("s", ss)
  fhist_surface = fhist_surface %>%
    as_tibble() %>%
    mutate(t = tt) %>% select(t, everything()) %>%
    gather(s, value, starts_with("s")) %>%
    arrange(t) %>%
    mutate(s = rep(ss, df$D),
           alpha = df$alpha_true,
           method = "fhist") %>%
    mutate(value = ifelse(s <= t, value, NA))
  
  
  load(f)
  # get surface for flode
  flode_surface = results$flode_results$surface[[2]]
  colnames(flode_surface) = paste0("s", ss)
  flode_surface = flode_surface %>%
    as_tibble() %>%
    mutate(t = tt) %>% select(t, everything()) %>%
    gather(s, value, starts_with("s")) %>%
    arrange(t) %>%
    mutate(s = rep(ss, df$D),
           alpha = df$alpha_true,
           method = "flode") %>%
  mutate(value = ifelse(s <= t, value, NA))
  

  bind_rows(true_surface, fhist_surface, flode_surface)
}


get_surfaces_genFhist = function(f){
  load(f)
  
  # get true surface
  df = as_tibble(results$params)
  tt = ss = seq(0, 1, length.out = df$D)
  true_surface = t(results$surface)
  colnames(true_surface) = paste0("s", ss)
  true_surface = true_surface %>% 
    as_tibble() %>%
    mutate(t = tt) %>% select(t, everything()) %>%
    gather(s, value, starts_with("s")) %>%
    arrange(t) %>%
    mutate(s = rep(ss, df$D),
           value = ifelse(s <= t, value, NA)) %>%  
    mutate(alpha = df$alpha_true,
           method = "truth")
  
  ## calculate surface for fhist
  fhist_surface = coef(results$fhist_results, n1 = df$D, n2 = df$D)$smterms$"ff(X1,s)"$coef %>%
    select(s = X1.smat, t = X1.tmat, value = value) %>%
    as_tibble()  %>%
    #mutate(s = round(s, 4)) #%>%
    mutate(#value = ifelse(s >= t, NA, value),
      alpha = df$alpha_true,
      method = "fhist") %>%
    select(s,t, value, alpha, method)
  
  ## interpolate surface to get on same grid as flode
  interp_obj = list(x = unique(fhist_surface$s), y = unique(fhist_surface$t),
                    z = matrix(fhist_surface$value, df$D, df$D, byrow = TRUE))
  
  loc_mat = list(x = ss, y = tt)
  fhist_surface = fields::interp.surface.grid(interp_obj, loc_mat)$z
  fhist_surface[df$D,] = fhist_surface[df$D-1,] # carry last value forward so not NA
  fhist_surface = fhist_surface * lower.tri(fhist_surface, diag = TRUE)# zero out upper triangle of the matrix
  
  colnames(fhist_surface) = paste0("s", ss)
  fhist_surface = fhist_surface %>%
    as_tibble() %>%
    mutate(t = tt) %>% select(t, everything()) %>%
    gather(s, value, starts_with("s")) %>%
    arrange(t) %>%
    mutate(s = rep(ss, df$D),
           alpha = df$alpha_true,
           method = "fhist") %>%
    mutate(value = ifelse(s <= t, value, NA))
  

  # get surface for flode
  flode_surface = results$flode_results$surface[[2]]
  colnames(flode_surface) = paste0("s", ss)
  flode_surface = flode_surface %>%
    as_tibble() %>%
    mutate(t = tt) %>% select(t, everything()) %>%
    gather(s, value, starts_with("s")) %>%
    arrange(t) %>%
    mutate(s = rep(ss, df$D),
           alpha = df$alpha_true,
           method = "flode") %>%
    mutate(value = ifelse(s <= t, value, NA))
  
  
  bind_rows(true_surface, fhist_surface, flode_surface)
}


get_predictions = function(f, f_fhist, new_data){
  load(f)
  params = as_tibble(results$params)

  alpha_vec = c(0.1, 0.5, 1, 4, 8, 12)
  if(params$alpha_true %in% alpha_vec){
    new_data = new_data[[which(params$alpha_true == alpha_vec)]]
  }else{
    tib = tibble(
      iter = params$iter,
      alpha_true = params$alpha_true,
      flode_mae = NA,
      flode_mspe = NA,
      fconc_mae = NA,
      fconc_mspe = NA,
      fhist_mae = NA,
      fhist_mspe = NA
    )
    return(tib)
  }
  
  # flode predictions
  flode = predict_flode(results$flode_results, 
                  new_data, 
                  forcing_functions = c("int", "x"))
  
  # fconc predictions
  xB_fconc = xB_fhist = matrix(NA, length(unique(new_data$trial)), 50)
  grid = sort(unique(new_data$time))
  trials = unique(new_data$trial)
  for(i in 1:length(trials)){
    xB = filter(new_data, trial == trials[i])$x * results$fconc_results$beta.pm$x
    xB_fconc[i,] = results$fconc_results$beta.pm$int + xB
  }
  
  # fhist predictions
  load(f_fhist)

  fhist_surf = coef(results$fhist_results, 
                    n1 = 50, n2 = 50)$smterms$"ff(x,s)"$coef %>%
    select(s = x.smat, t = x.tmat, value = value) %>% as_tibble() 
  
  surface_x = make_surface_fhist(fhist_surf, label = "x_fhist")
  b0 = coef(results$fhist_results)$pterms[1]
  beta0_fhist = as.numeric(coef(results$fhist_results, 
                                n1 = 50, n2 = 50)$smterms$"Intercept(t)"$coef$value)
  
  for(i in 1:length(trials)){
    x = filter(new_data, trial == trials[i])$x
    xb = integrate_surface(grid, surface_x, covar = x) 
    
    xB_fhist[i,] = b0 + beta0_fhist + xb
  }
  
  # aggregate results
  message(params$iter)
  
  tibble(
    iter = params$iter,
    alpha_true = params$alpha_true,
    flode_mae = mean(abs(flode$value -flode$yhat_flode)),
    flode_mspe = mean((flode$value -flode$yhat_flode)^2),
    fconc_mae = mean(abs(new_data$value - as.vector(t(xB_fconc)))),
    fconc_mspe = mean((new_data$value - as.vector(t(xB_fconc)))^2),
    fhist_mae = mean(abs(new_data$value - as.vector(t(xB_fhist)))),
    fhist_mspe = mean((new_data$value - as.vector(t(xB_fhist)))^2)
  )
}
  



get_predictions_genKalman = function(f, new_data){
  load(f)
  params = as_tibble(results$params)
  
  dat = new_data$data
  dat = dat %>% mutate(int = 1) %>%
    group_by(trial) %>%
    mutate(y0 = first(value))  %>%
    ungroup() 


  # flode predictions
  flode = predict_flode(results$flode_results, 
                        dat, 
                        forcing_functions = c("int", "x"))  
  
  # fhist predictions
  xB_fhist = matrix(NA, length(unique(dat$trial)), 50)
  grid = sort(unique(dat$time))
  trials = unique(dat$trial)
  
  fhist_surf = coef(results$fhist_results, 
                    n1 = 50, n2 = 50)$smterms$"ff(x,s)"$coef %>%
    select(s = x.smat, t = x.tmat, value = value) %>% as_tibble() 
  
  surface_x = make_surface_fhist(fhist_surf, label = "x_fhist")
  b0 = coef(results$fhist_results)$pterms[1]
  beta0_fhist = as.numeric(coef(results$fhist_results, 
                                n1 = 50, n2 = 50)$smterms$"Intercept(t)"$coef$value)
  
  for(i in 1:length(trials)){
    x = filter(dat, trial == trials[i])$x
    xb = integrate_surface(grid, surface_x, covar = x) 
    
    xB_fhist[i,] = b0 + beta0_fhist + xb
  }
  
  # aggregate results
  message(params$iter)
  
  tibble(
    iter = params$iter,
    flode_mae = mean(abs(flode$value -flode$yhat_flode)),
    flode_mspe = mean((flode$value -flode$yhat_flode)^2),
    fhist_mae = mean(abs(dat$value - as.vector(t(xB_fhist)))),
    fhist_mspe = mean((dat$value - as.vector(t(xB_fhist)))^2)
  )
}
  



get_predictions_genFhist = function(f, new_data){
  load(f)
  params = as_tibble(results$params)
  
  dat = new_data$data
  df = new_data$df %>% mutate(int = 1) %>%
    group_by(trial) %>%
    mutate(y0 = first(value))  %>%
    ungroup() 
  
  # flode predictions
  flode = predict_flode(results$flode_results, 
                        df, 
                        forcing_functions = c("int", "x"))
  
  # fconc predictions
  xB_fconc = xB_fhist = matrix(NA, length(unique(df$trial)), 50)
  grid = sort(unique(df$time))
  trials = unique(df$trial)
  for(i in 1:length(trials)){
    xB = filter(df, trial == trials[i])$x * results$fconc_results$beta.pm$x
    xB_fconc[i,] = results$fconc_results$beta.pm$int + xB
  }
  
  
  
  # fhist predictions
  fhist_surf = coef(results$fhist_results, 
                    n1 = 50, n2 = 50)$smterms$"ff(X1,s)"$coef %>%
    select(s = X1.smat, t = X1.tmat, value = value) %>% as_tibble() 
  
  surface_x = make_surface_fhist(fhist_surf, label = "x_fhist")
  b0 = coef(results$fhist_results)$pterms[1]
  beta0_fhist = as.numeric(coef(results$fhist_results, 
                                n1 = 50, n2 = 50)$smterms$"Intercept(t)"$coef$value)
  
  
  for(i in 1:length(trials)){
    x = filter(df, trial == trials[i])$x
    xb = integrate_surface(grid, surface_x, covar = x) 
    
    xB_fhist[i,] = b0 + beta0_fhist + xb
  }
  
  # aggregate results
  message(params$iter)
  
  tibble(
    iter = params$iter,
    alpha_true = params$alpha_true,
    flode_mae = mean(abs(flode$value -flode$yhat_flode)),
    flode_mspe = mean((flode$value -flode$yhat_flode)^2),
    fconc_mae = mean(abs(df$value - as.vector(t(xB_fconc)))),
    fconc_mspe = mean((df$value - as.vector(t(xB_fconc)))^2),
    fhist_mae = mean(abs(df$value - as.vector(t(xB_fhist)))),
    fhist_mspe = mean((df$value - as.vector(t(xB_fhist)))^2)
  )
}


