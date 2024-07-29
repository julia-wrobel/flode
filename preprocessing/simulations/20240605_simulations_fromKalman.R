####################################################################
# Julia Wrobel
# June 5, 2024
# 
# This file produces simulations for flode paper for generating data from kalman model (supp Figure 6)
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

source(here::here("preprocessing", "simulated_paw_kalman.R"))
source(here::here("preprocessing", "estimate_flode.R"))
source(here::here("preprocessing", "make_pffr.R"))
source(here::here("preprocessing", "utils.R"))
source(here::here("preprocessing", "search_alpha.R"))

###############################################################
## set simulation design elements 
###############################################################
N_vec = c(100)
D_vec = c(50)
a = c(0.001)

seed_start = 1000 
N_iter = 50
miter = 100


params = expand.grid(seed_start = seed_start, 
                     N_vec = N_vec,
                     D_vec = D_vec)


## record date for analysis; create directory for results
Date = gsub("-", "", Sys.Date())


###############################################################
## start simulation code
###############################################################


for(scenario in 1:dim(params)[1]){
  
  results_mat = matrix(NA, nrow = 1, ncol = 14)
  colnames(results_mat) = c("iter", "seed", "alpha_est", "sse_flode", 
                            "maxiter", "N", "D","surfaceErr_flode", "time_flode", 
                            "sse_kalman", "surfaceErr_kalman", "time_kalman", "sigma_flode", "sigma_kalman")
  
  
  ###############################################################
  ## set simulation design elements 
  ###############################################################
  N = params$N_vec[scenario]
  D = params$D_vec[scenario]
  SEED.START = params$seed_start[scenario]
  
  
  results_mat[, 6] = N
  results_mat[, 7] = D
  
  results = foreach(iter = 1:N_iter) %dopar% {
    
    seed.iter = (SEED.START - 1)*N_iter + iter
    results_mat[, 3] = seed.iter
    
    ## generate data
    simulated_data = simulate_kalman(I = N, 
                                    D = D, 
                                    sigma = 0.1, 
                                    rand_int = TRUE, 
                                    lambda0 = 10,
                                    seed = seed.iter)
    
    
    ################################################################################
    ## flode 
    dat = simulated_data$data %>% mutate(int = 1) %>%
      group_by(trial) %>%
      mutate(y0 = first(value))  %>%
      ungroup() 
    
    y0_df = simulated_data$data %>% 
      group_by(trial) %>%
      slice(1) %>%
      ungroup() 
    
    set.seed(seed.iter)
    
    tic(quiet = TRUE)
    
    # grid search to initialize alpha
    alpha_start = search_alpha(data = dat, Kt = 20)$best_alpha
    
    # run flode
    flode_results = estimate_flode(dat, 
                                   Kt = 20, 
                                   alpha0 = alpha_start,  
                                   max_iter = miter, 
                                   forcing_functions = c("int", "x"),
                                   tol =  0.0001, 
                                   sigma_b = 10,
                                   sigma_d = 10,
                                   alpha_upper = 20,
                                   random_effects = TRUE,
                                   initial_position = TRUE,
                                   y0_vec = y0_df$value)
    time_flode = toc()
    
    results_mat[, 1] = iter
    
    results_mat[, 3] = flode_results$alpha
    results_mat[, 4] = flode_results$sse
    
    ################################################################################
    ## get surface errors and computation time for flode
    
    # flode surface error
    results_mat[, 8] = sum((simulated_data$surface - flode_results$surface[[2]])^2)
    
    # get computation time
    results_mat[, 9] = time_flode$toc -time_flode$tic
    results_mat[, 13] = flode_results$sigma
    
    colnames(results_mat) = c("iter", "seed", "alpha_est", "sse_flode", 
                              "maxiter", "N", "D","surfaceErr_flode", 
                              "time_flode", "sse_kalman", "surfaceErr_kalman", "time_kalman", 
                              "sigma_flode", "sigma_kalman")
    
    ################################################################################
    ## kalman
    tic(quiet = TRUE)
    kalman_results_re = make_pffr(dat, random_int = TRUE)
    time_kalman = toc()
    

    results = list(params = results_mat,
                   times = list(time_flode, time_kalman),
                   flode_results = flode_results,
                   fhist_results = kalman_results_re)
    
    
    filename = paste0("/Users/JWROBE8/onedrive/Data/flode/kalman/", Date, "_iter_", iter, "_scenario_" ,scenario, ".RDA")
    save(results,
         file = filename)
    
  } # end N_iter foreach loop
  
  message(paste0("parameter scenario: ", scenario))

} # end scenario loop 

###############################################################
## end sim
###############################################################


