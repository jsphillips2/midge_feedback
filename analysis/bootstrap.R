options(mc.cores = parallel::detectCores()-2)

bootMer0 <- function(mod, term, nAGQ=0, nboot=1, verbose=F, c=1) {
  
  formula <- as.formula(paste0(".~.-",term))
  
  mod <- update(mod, control = glmerControl(calc.derivs=F))
  mod0 <- update(mod, formula=formula)
  
  z.true <- logLik(mod) - logLik(mod0)
  sim.data <- model.frame(mod)
  z<- parallel::mclapply(1:nboot, function(j){
    simY <- simulate(mod0)
    if (c == 1) {
      sim.data$Y1 <- simY[,1]
    }
    if (c == 2) {
      sim.data$Y1 <- simY[,1][,1]
      sim.data$Y2 <- simY[,1][,2]
    }
    mod.boot <- update(mod, data=sim.data, start=list(theta=attributes(mod)$theta))
    mod0.boot <- update(mod0, data=sim.data, start=list(theta=attributes(mod0)$theta))
    
    z <- logLik(mod.boot) - logLik(mod0.boot)
    return(z)
  }) %>% unlist()
  P <- mean(z.true < z)
  Chi <- 2*z.true
  
  return(tibble(term = term, p = P, chi = 2*z.true))
}

bootLRT <- function(mod, nAGQ=0, nboot=1, verbose=F, c=1) {
  terms <- attributes(terms(mod))$term.labels
  d <- lapply(terms, function(x){
    bootMer0(mod = mod, term = x, nAGQ = nAGQ, nboot = nboot, verbose = verbose, c = c)}) %>% 
    bind_rows()
  return(d)
}





