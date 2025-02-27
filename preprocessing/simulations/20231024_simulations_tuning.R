####################################################################
# Julia Wrobel
# October 24, 2023
# 
# This file produces simulations for flode paper, specifically for varying over tuning parameter values.
# Compares flode to functional historical regression
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
source(here::here("preprocessing", "utils.R"))
source(here::here("preprocessing", "search_alpha.R"))

###############################################################
## set simulation design elements 
###############################################################
N_vec = c(100)
D_vec = c(50)
lambda_d_vec = c(1, 10, 100)
alpha_vec = c(4)
a = c(0.01, 0.001)

seed_start = 1000 
Kt = c(10, 20, 40)
N_iter = 25
miter = 200


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
	
	results_mat = matrix(NA, nrow = 1, ncol = 16)
	colnames(results_mat) = c("iter","scenario", "seed", "alpha_est", "lambda_d_est", "sse_flode", "logLik", 
													"a", "N", "D", "alpha_true", "lambda_d_true", "surfaceErr_flode", "time_flode", 
													"Kt", "lambda_b")

	results_mat[, 2] = scenario

	
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
	
	results_mat[, 8] = a
	results_mat[, 9] = N
	results_mat[, 10] = D
	results_mat[, 11] = alpha0
	results_mat[, 12] = lambda_d0
	results_mat[, 15] = Kt
	results_mat[, 16] = lambda_d0
	
	results = foreach(iter = 1:N_iter, .combine = 'rbind') %dopar% {
		
		seed.iter = (SEED.START - 1)*N_iter + iter
		results_mat[, 3] = seed.iter
		
		## generate data
		simulated_data = simulate_flode(I = N, 
																		D = D, 
																		sigma_y0 = 5, 
																		sigma = 0.1, 
																		alpha = alpha0, 
																		rand_int = TRUE, 
																		lambda0 = 20,
																		seed = seed.iter)
		
		
		################################################################################
		## flode 
		dat = simulated_data$data %>% mutate(int = 1)
		
		set.seed(seed.iter)
		
		tic(quiet = TRUE)
		
		# grid search to initialize alpha
		alpha_start = search_alpha(data = dat, Kt = Kt)$best_alpha
		
		# run flode
		flode_results = estimate_flode(dat, 
																	 Kt = Kt, 
																	 alpha0 = alpha_start,  
																	 max_iter = miter, 
																	 forcing_functions = c("int", "x"),
																	 tol =  0.00001, 
																	 lambda_b = lambda_d0,
																	 lambda_d = lambda_d0,
																	 random_effects = TRUE,
																	 initial_position = TRUE,
																	 y0_vec = simulated_data$y0,
																	 a = a)
		time_flode = toc()
		
		results_mat[, 1] = iter
		
		results_mat[, 4] = flode_results$alpha
		results_mat[, 5] = flode_results$lambda_d
		results_mat[, 6] = flode_results$sse
		results_mat[, 7] = flode_results$logLikelihood
		
		################################################################################
		## get surface errors and computation time for flode
		
		# flode surface error
		results_mat[, 13] = sum((simulated_data$surface - flode_results$surface[[2]])^2)
		
		# get computation time
		results_mat[, 14] = time_flode$toc -time_flode$tic
		
		
		results_mat
	} # end N_iter foreach loop

	message(paste0("parameter scenario: ", scenario))
	
	filename = paste0(here::here("output", "simulation_results"), "flode_tuning_", scenario, ".RDA")
	save(results,
			 file = filename)
} # end scenario loop 

###############################################################
## end sim
###############################################################


