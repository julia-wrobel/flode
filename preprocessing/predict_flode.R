predict_flode = function(flode_model, 
                         new_data,
                         forcing_functions = c(""),
                         boot = FALSE){
  
  if(boot){
    new_data = as_tibble(new_data) %>% 
      mutate(old_trial = trial,
             trial = row_number()) %>%
      unnest(data) 
  }
  
  x_data = new_data %>% select(trial, all_of(forcing_functions))  %>% nest(data = c(-trial)) %>% pull(data)
  basis_time = filter(new_data, trial == first(trial))$time # assumes grid is even across subjects
  spline_basis = flode_model$spline_basis
  
  N = length(unique(new_data$trial))
  D = length(basis_time)
  alpha_cur = flode_model$alpha
  sigma_cur = flode_model$sigma
  beta_coefs = flode_model$beta_coefs
  P = dim(beta_coefs)[2]
  Kt = dim(beta_coefs)[1]
  Pen = flode_model$var_components$Pen
  
  # make y* term
  y0_vec = (new_data %>% group_by(trial) %>% slice(1) %>% ungroup())$value # data should be sorted by time for each subject
  y0_cur = rep(y0_vec, each = D)
  new_data = new_data %>% mutate(Yo = y0_cur, 
                                 y0_star = Yo * exp(- alpha_cur * time))
  
  thetaStar = integrate_alpha(basis_time, spline_basis, alpha_cur, covar = NULL)
  xstar_ls = lapply(x_data, make_Xstar, alpha = alpha_cur, spline_basis = spline_basis, basis_time = basis_time,
                    thetaStar = thetaStar)
  xstar = do.call(rbind, xstar_ls)
  new_data = new_data %>% mutate(yhat_flode = y0_star + as.numeric(xstar %*% as.vector(beta_coefs)))
  
  if(boot){
    Dstar = thetaStar
    sigma_b_cur = flode_model$sigma_b
    sigma_d_cur = flode_model$sigma_d
 
    iter = 0
    trial_ids = unique(new_data$trial)
    delta_coefs = matrix(0, nrow = Kt, ncol = N)
    deltaP_squared = matrix(0, nrow = 1, ncol = N)
    beta_coefs = solve(crossprod(xstar) + sigma_cur / (2 * sigma_b_cur) * Matrix::bdiag(rep(list(Pen + t((Pen))), P))) %*% (t(xstar) %*% (new_data$value - new_data$y0_star - as.vector(Dstar %*% delta_coefs)))
    beta_coefs = matrix(beta_coefs, nrow = Kt, ncol = P)
    
    while(iter < 20){
      C = solve(1/sigma_d_cur * Pen + crossprod(Dstar)/sigma_cur) 
      trace_PC = sum(diag(Pen %*% C))
      C_Dstar = C %*% t(Dstar)
      
      for(i in 1:length(trial_ids)){
        trial_df = filter(new_data, trial == trial_ids[i])
        
        delta_coefs[, i] =   (C_Dstar %*% (trial_df$value - trial_df$y0_star -
                                             xstar_ls[[i]] %*% as.vector(beta_coefs)))/sigma_cur
        deltaP_squared[, i] = trace_PC + crossprod(delta_coefs[,i], Pen %*% delta_coefs[,i])
      }
      
      # update sigma
      new_data = new_data %>% mutate(delta =  as.vector(spline_basis %*% delta_coefs),
                                     deltaStar = as.vector(Dstar %*% delta_coefs),
                                     yhat = y0_star + as.numeric(xstar %*% as.vector(beta_coefs)) + deltaStar,
                                     err = value - yhat)
      RSS = sum(new_data$err^2)
      sigma_cur = as.numeric(RSS)/(N * D)
      
      # update sigma_d
      sigma_d_cur = sum(deltaP_squared)/(N * Kt)
      # update sigma_b
      sigma_b_cur = sum(diag(crossprod(beta_coefs, Pen %*% beta_coefs)))/(Kt*P)
      
      #
      beta_coefs = solve(crossprod(xstar) + sigma_cur / (2 * sigma_b_cur) * Matrix::bdiag(rep(list(Pen + t((Pen))), P))) %*% (t(xstar) %*% (new_data$value - new_data$y0_star - as.vector(Dstar %*% delta_coefs)))
      beta_coefs = matrix(beta_coefs, nrow = Kt, ncol = P)

      iter = iter + 1
    }
    
    
    beta_cur = spline_basis %*% beta_coefs
    colnames(beta_cur) = paste0("beta", 0:(P-1))
    beta_cur = as_tibble(beta_cur) %>% mutate(time = basis_time) 
        
    beta_cur
  }else{
    new_data
  }
  
  
}



# fhist surfaces
make_surface_fhist = function(fhist_surf, D = 50, label = ""){
  ## interpolate surface to get on same grid as flode
  interp_obj = list(x = unique(fhist_surf$s), y = unique(fhist_surf$t),
                    z = matrix(fhist_surf$value, D, D, byrow = TRUE))
  
  t = s = seq(0, 1, length.out = D)
  loc_mat = list(x = s, y = t)
  surface = fields::interp.surface.grid(interp_obj, loc_mat)$z
  surface[D,] = surface[D-1,] # carry last value forward so not NA
  surface =	surface * lower.tri(surface, diag = TRUE) %>% as_tibble() 
  
  
  surface %>% 
    mutate(t = t) %>%
    pivot_longer(V1:V50, names_to = "s", values_to = "value") %>%
    mutate(label = label,
           s = rep(seq(0, 1, length.out = D), 50)) 
  
}



integrate_surface = function(basis_time, surface, covar){
  D = length(basis_time)
  
  beta_covar_f = function(s0, t0 = 1, basis_time, surface, covar){
    filter(surface, s == s0, t == t0)$value * covar[which(basis_time == s0)] 
  }
  
  Gstar = sapply(basis_time, FUN = function(t){cumtrapz(basis_time, beta_covar_f(basis_time, t, basis_time, surface,covar))[D,]})
  return(t(Gstar))
}



# integrate_surface = function(basis_time, surface, covar){
#   D = length(basis_time)
#   
#   beta_covar_f = function(s0, t0 = 1, basis_time, surface, covar){
#     filter(surface, s < t0, t == t0)$value * covar[which(basis_time < t0)] 
#   }
#   
#   Gstar = sapply(basis_time, FUN = function(t){cumtrapz(basis_time, beta_covar_f(basis_time, t, basis_time, surface,covar))[D,]})
#   return(t(Gstar))
# }
