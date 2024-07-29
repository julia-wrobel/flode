pffr_sim = function(scenario = c("int", "ff"), n = 100, nxgrid = 50,
          SNR = 10,  limits = TRUE, 
          sigma_d = 10,
          rand_int = TRUE){
  mc <- match.call()
  for (i in 2:length(mc)) if (is.symbol(mc[[i]])) 
    mc[[i]] <- get(deparse(mc[[i]]), envir = parent.frame())
  rf <- function(x = seq(0, 1, length = 100), bs.dim = 7, center = FALSE) {
    nk <- bs.dim - 2
    xu <- max(x)
    xl <- min(x)
    xr <- xu - xl
    xl <- xl - xr * 0.001
    xu <- xu + xr * 0.001
    dx <- (xu - xl)/(nk - 1)
    kn <- seq(xl - dx * 3, xu + dx * 3, length = nk + 4 + 
                2)
    X <- splines::spline.des(kn, x, 4, x * 0)$design
    drop(X %*% rnorm(bs.dim))
  }
  s <- seq(0, 1, length = nxgrid)
  t <- seq(0, 1, length = nxgrid)
  mu.t <- matrix(1 + dbeta(t, 2, 7), nrow = n, ncol = nxgrid, 
                 byrow = TRUE)
  data <- list()
  etaTerms <- list()
  etaTerms$int <- mu.t
  data$X1 <- I(t(replicate(n, rf(s))))
  L <- matrix(1/nxgrid, ncol = nxgrid, nrow = n)
  LX1 <- L * data$X1
  ## make coefficient surface
  surface_df =
    expand_grid(
      s = s,
      t = t
    ) |> 
    mutate(
      coef = map2_dbl(s, t, \(x,y) dmvnorm(x = c(x, y), mean = c(.25, .75), sigma = matrix(c(.0075, 0, 0, .0075), 2, 2))),
    )
  
  if(limits){
    surface_df = surface_df %>% mutate(coef = ifelse(s < t, coef, 0))
  }
  
  beta1.st = matrix(surface_df$coef, ncol = nxgrid, 
                    nrow = nxgrid, byrow = TRUE)
  #beta1.st <- outer(s, t, test1)

  etaTerms$X1 <- LX1 %*% beta1.st
  data$xlin <- I(rnorm(n))
  beta.t <- matrix(scale(-dnorm(4 * (t - 0.2))), nrow = n, 
                   ncol = nxgrid, byrow = T)
  etaTerms$xlin <- data$xlin * beta.t
  data$xsmoo <- I(rnorm(n))
  etaTerms$xsmoo <- outer(drop(scale(cos(data$xsmoo))), (t - 
                                                           0.5), "*")
  data$xfactor <- sample(gl(3, n/3), replace = TRUE)
  etaTerms$xfactor <- 2 * as.numeric(data$xfactor) + sin(2 * 
                                                           outer(as.numeric(data$xfactor), t))
  if ("2factor" %in% scenario) {
    data$x2factor <- sample(gl(3, n/3), replace = TRUE)
    etaTerms$x2factor <- 2 * as.numeric(data$x2factor) + 
      cos(2 * outer(as.numeric(data$x2factor), t))
  }

  if (length(scenario) == 1) {
    eta <- mu.t + switch(scenario, int = 0, all = Reduce("+", 
                                                         etaTerms), ff = etaTerms$X1, 
                         lin = etaTerms$xlin, 
                         smoo = etaTerms$xsmoo, te = etaTerms$xte, const = etaTerms$xconst, 
                         factor = etaTerms$xfactor)
  }
  else {
    stopifnot(all(scenario %in% c("int", "ff", "lin", "smoo", 
                                  "te", "const", "factor", "2factor")))
    eta <- 0 * mu.t
    if ("int" %in% scenario) 
      eta <- eta + mu.t
    if ("ff" %in% scenario) 
      eta <- eta + etaTerms$X1
    if ("lin" %in% scenario) 
      eta <- eta + etaTerms$xlin
    if ("smoo" %in% scenario) 
      eta <- eta + etaTerms$xsmoo
    if ("te" %in% scenario) 
      eta <- eta + etaTerms$xte
    if ("const" %in% scenario) 
      eta <- eta + etaTerms$xconst
    if ("factor" %in% scenario) 
      eta <- eta + etaTerms$xfactor
    if ("2factor" %in% scenario) 
      eta <- eta + etaTerms$x2factor
  }
  eps <- sd(as.vector(eta))/sqrt(SNR) * matrix(scale(rnorm(n * 
                                                             nxgrid)), nrow = n)
  if(rand_int){
    sigma_d = rep(sigma_d, 10 * n)
    beta_basis = bs(t, df = 10, intercept = TRUE)
    re_coef = matrix(rnorm(n * 10, 0, sqrt(sigma_d)), nrow = n, ncol = 10)
    re_mat = re_coef %*% t(beta_basis) 
  }else{
    re_mat = 0
  }
  
  Y <- eta + eps + re_mat
  
  
  data = as.data.frame(data, 
                       rownames = 1:n)
  data$Y = I(Y)
  data$trial = factor(1:n)
  
  
  df = tibble(trial = rep(data$trial, each = nxgrid),
               time = rep(t, n),
               value = as.vector(t(Y)),
               x = as.vector(t(data$X1)),
               xB = as.vector(t(etaTerms$X1)),
               re = as.vector(t(re_mat)),
  )
  return(list(data = data, 
              df = df,
              s = s, 
              t = t, 
              surface = beta1.st,
              surface_df = surface_df,
              truth = list(eta = eta, 
                           etaTerms = etaTerms)))
}
