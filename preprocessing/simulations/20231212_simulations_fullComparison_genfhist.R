####################################################################
# Julia Wrobel
# 
# This file produces simulations for flode paper.
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
library(vbvs.concurrent)
library(refund)


# get parallelization set up
library(foreach)
library(doParallel)
n_cores = detectCores() - 1
registerDoParallel(n_cores)

###############################################################
## define or source functions used in code below
###############################################################

source(here::here("preprocessing", "pffr_sim.R"))
source(here::here("preprocessing", "estimate_flode.R"))
source(here::here("preprocessing", "utils.R"))
source(here::here("preprocessing", "search_alpha.R"))



###############################################################
## set simulation design elements 
###############################################################

D_vec = c(50)
lambda_d_vec = c(10)
alpha_vec = c(12)
a = c(0.001)

seed_start = 1000 
Kt = c(10)

N_vec = c(100)
N_iter = 50 # 50 
miter = 200 #100



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
	
	
	results_mat = matrix(NA, nrow = 1, ncol = 5)
	colnames(results_mat) = c("iter", "seed", "N", "D", "alpha_true")
	
	results_mat[, 3] = N
	results_mat[, 4] = D
	results_mat[, 5] = alpha0
	
	results = foreach(iter = 1:N_iter) %dopar% {
		
		seed.iter = (SEED.START - 1)*N_iter + iter
		# set iteration 
		results_mat[, 1] = iter
		results_mat[, 2] = seed.iter
		
		## generate data
		set.seed(seed.iter)
		simulated_data = pffr_sim(scenario = c("int", "ff"),
		                          N, nxgrid = D, SNR = 10,
		                          sigma_d = 1,
		                          rand_int = TRUE)
		
		
		dat = simulated_data$data
		df = simulated_data$df %>% mutate(int = 1) %>%
		  group_by(trial) %>%
		  mutate(y0 = first(value))  %>%
		  ungroup() 
		
		y0_df = df %>% 
		  group_by(trial) %>%
		  slice(1) %>%
		  ungroup() 
		
		################################################################################
		## flode 
		################################################################################

		tic(quiet = TRUE)
		alpha_start = search_alpha(data = df,
		                           Kt = Kt,
		                           alpha_min = 4,
		                           alpha_max = 30)$best_alpha
		flode_results = estimate_flode(df,
		                               Kt = Kt,
		                               alpha0 = alpha_start,
		                               max_iter = miter,
		                               forcing_functions = c("int", "x"),
		                               tol =  0.0001,
		                               sigma_b = lambda_d0,
		                               sigma_d = lambda_d0,
		                               random_effects = TRUE,
		                               initial_position = TRUE,
		                               y0_vec = simulated_data$y0,
		                               alpha_upper = alpha_start + 5)
		time_flode = toc()

		################################################################################
		## fhist
		################################################################################
		tic(quiet = TRUE)
		s = t = simulated_data$t
		fhist_results = pffr(Y ~ s(trial, bs = "re") + ff(X1, xind = s,limits = "s<t",
		                                        splinepars = list(bs = "ps", 
		                                                          m = list(c(2, 1), c(2, 1)), 
		                                                          k= c(15, 15))), 
		           yind = t, dat)
		time_fhist = toc()
		
		################################################################################
		## fconc
		################################################################################

		tic(quiet = TRUE)
		fconc_results = vb_concurrent(value ~ int + x | time,
		                              id.var = "trial",
		                              data = df,
		                              Kt = Kt,
		                              Kp = 15,
		                              t.min = 0, t.max = 1,
		                              standardized = TRUE)
		time_fconc = toc()
		################################################################################
		## end modeling
		################################################################################
		
		
	 ## save folds
		
		results = list(params = results_mat,
		     times = list(time_flode, time_fhist, time_fconc),
		     fhist_results = fhist_results,
		     fconc_results = fconc_results,
		     flode_results = flode_results,
		     surface = simulated_data$surface)
		
		filename = paste0("/Users/JWROBE8/Desktop/flode/", Date, "_iter_", iter, "_scenario_" ,scenario, ".RDA")
		#filename = paste0("/Users/JWROBE8/onedrive/Data/flode/funreg_comparison/genfhist/", Date, "_iter_", iter, "_scenario_" ,scenario, ".RDA")
		save(results,
		     file = filename)
	} # end N_iter foreach loop

	message(paste0("parameter scenario: ", scenario))
	
	#filename = paste0(here::here("output", "simulation_results", "funreg_comparison"), "/funregComparison_", Date, "_",scenario, ".RDA")
	#save(results,
	#		 file = filename)
} # end scenario loop 

###############################################################
## end sim
###############################################################


