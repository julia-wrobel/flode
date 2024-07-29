make_pffr_pUp = function(data, covar = "x", random_int = TRUE){
	N = length(unique(data$trial))
	D = length(filter(data, trial == first(trial))$time)
	
	y_fhist = matrix(data$value, nrow = N, ncol = D, byrow = TRUE)
	x_fhist1 = matrix(data[[covar[1]]], nrow = N, ncol = D, byrow = TRUE)
	x_fhist2 = matrix(data[[covar[2]]], nrow = N, ncol = D, byrow = TRUE)
	
	y0 = data %>%
		group_by(trial) %>%
		summarize(y0 = first(y0)) %>% ungroup() %>% pull(y0)
	
	fhist_df = data.frame(trial = factor(1:N), y0 =  y0)
	fhist_df$y = y_fhist
	fhist_df$x1 = x_fhist1
	fhist_df$x2 = x_fhist2
	t = s = seq(0, 1, length.out = D)
	
	if(random_int){
		mod = pffr(y ~ s(trial, bs = "re", k = 25) + ff(x1, xind = s, limits = "s<t",
		                                                splinepars = list(bs = "ps", 
		                                                                  m = list(c(2, 1), c(2, 1)), 
		                                                                  k= c(15, 15))) + 
		             ff(x2, xind = s, limits = "s<t",
		                splinepars = list(bs = "ps", 
		                                  m = list(c(2, 1), c(2, 1)), 
		                                  k= c(15, 15))), 
							 yind = t, fhist_df)
	}else{
		mod = pffr(y ~ y0 + ff(x, xind = s, limits = "s<t"), yind = t, fhist_df)
	}
	return(list(mod = mod, fhist_df = fhist_df)) 
}