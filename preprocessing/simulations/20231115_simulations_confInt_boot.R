####################################################################
# Julia Wrobel
# 
# This file produces simulations for flode paper for getting coverage for bootstrap confidence intervals
# 
####################################################################

library(splines)
library(tibble)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(mvtnorm)
library(pracma)
library(tictoc)
library(refund)
library(tidymodels)


# get parallelization set up
library(foreach)
library(doParallel)
n_cores = detectCores() - 1
registerDoParallel(n_cores)

###############################################################
## define or source functions used in code below
###############################################################

source(here::here("preprocessing", "simulated_paw.R"))
source(here::here("preprocessing", "estimate_flode.R"))
source(here::here("preprocessing", "predict_flode.R"))
source(here::here("preprocessing", "utils.R"))
source(here::here("preprocessing", "search_alpha.R"))

###############################################################
## set simulation design elements 
###############################################################
N_vec = c(100)
D_vec = c(50)
lambda_d_vec = c(10)
alpha_vec = c(4)
a = c(0.001)

seed_start = 1000 
Kt = c(10)

N_iter = 500
miter = 100
n_boot = 500


params = expand.grid(seed_start = seed_start, 
										 N_vec = N_vec,
										 D_vec = D_vec, 
										 lambda_d_vec = lambda_d_vec,
										 alpha_vec = alpha_vec,
										 a = a,
										 Kt = Kt)


## record date for analysis; create directory for results
Date = gsub("-", "", Sys.Date())


###############################################################
## start simulation code
###############################################################


## need to parallelize this part
## actually, parallelize the iterations?
for(scenario in 1:dim(params)[1]){
	

	
	###############################################################
	## set simulation design elements 
	###############################################################
	a = params$a[scenario]
	N = params$N_vec[scenario]
	D = params$D_vec[scenario]
	SEED.START = params$seed_start[scenario]
	alpha0 = params$alpha_vec[scenario]
	lambda_d0 = params$lambda_d_vec[scenario]
	Kt = params$Kt[scenario]
	
	
	results_mat = matrix(NA, nrow = 1, ncol = 6)
	colnames(results_mat) = c("iter", "seed", "N", "D", "alpha_true", "time_flode")
	
	results_mat[, 3] = N
	results_mat[, 4] = D
	results_mat[, 5] = alpha0
	
	results = foreach(iter = 1:N_iter) %dopar% {
		
		seed.iter = (SEED.START - 1)*N_iter + iter
		results_mat[, 2] = seed.iter
		
		## generate data
		simulated_data = simulate_flode(I = N, 
																		D = D, 
																		sigma_y0 = 5, 
																		sigma = 0.05, 
																		alpha = alpha0, 
																		rand_int = TRUE, 
																		lambda0 = 20,
																		autoregressive = TRUE,
																		seed = seed.iter)
		
		
		################################################################################
		## flode 
		dat = simulated_data$data %>% mutate(int = 1)
		
		set.seed(seed.iter)
		
		tic(quiet = TRUE)
		
		# run flode
		flode_results = estimate_flode(dat, 
																	 Kt = Kt, 
																	 alpha0 = alpha0,  
																	 max_iter = miter, 
																	 forcing_functions = c("int", "x"),
																	 tol =  0.001, 
																	 sigma_b = lambda_d0,
																	 sigma_d = lambda_d0,
																	 random_effects = TRUE,
																	 initial_position = TRUE,
																	 y0_vec = simulated_data$y0,
																	 a = a)
		time_flode = toc()
		
		results_mat[, 1] = iter
		# get computation time
		results_mat[, 6] = time_flode$toc -time_flode$tic

		# do bootstrap
		dat = dat %>%
		  select(trial, int, time, value, x, y0) %>%
		  nest(data = c(-trial))
		
		
		dat_boot <- bootstraps(dat,
		                       times = n_boot)
		
		
		dat_betas <- dat_boot %>%
		  mutate(beta = map(splits, ~predict_flode(flode_model = flode_results,
		                                           new_data = .,
		                                           forcing_functions = c("int",  "x"),
		                                           boot = TRUE) )
		  )
		
		
		beta_full = flode_results$beta %>%
		  pivot_longer(beta0:beta1, names_to = "coef", values_to = "beta_hat")
		
		beta_true = simulated_data$coefficient_fns %>%
		  select(time, beta0, beta1) %>%
		  pivot_longer(beta0:beta1, names_to = "coef", values_to = "beta_true")
		
		ci_df = dat_betas %>%
		  unnest(beta) %>%
		  pivot_longer(beta0:beta1, names_to = "coef", values_to = "value") %>%
		  group_by(time, coef) %>%
		  summarize(lower = quantile(value, 0.025),
		            upper = quantile(value, 0.975),
		            sd = sd(value))  %>%
		  ungroup() %>%
		  mutate(beta_hat = beta_full$beta_hat,
		         beta_true = beta_true$beta_true,
		         lower_wald = beta_hat - 1.96 * sd,
		         upper_wald = beta_hat + 1.96 * sd,
		         iter = iter,
		         seed = seed.iter)
		
		message(paste0("iteration: ", iter))
		## save wald boot and quantile boot
		list(params = results_mat,
		     CI = ci_df)
		
	} # end N_iter foreach loop

	message(paste0("parameter scenario: ", scenario))
	
	filename = paste0(here::here("output", "simulation_results"), "/flode_confInt_bootAR1_", Date, "_",scenario, ".RDA")
	save(results,
			 file = filename)
} # end scenario loop 

###############################################################
## end sim
###############################################################


