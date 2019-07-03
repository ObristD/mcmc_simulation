# **************** Monte Carlo Markov Chain Simulation (MCMC) ***********************************
# @ Author: Dario Obrist
#
# @ Method: Metropolis Hastings Algorithm (randrom walk)
#
# @ Param:
## N:             number of simulations
## FUN:           function 
## ...:           additional patameters
## burnin:        BurnIn phase (number of simulations to be cut of)
## x0:            startvalue
## discrete:      bool for discrete distribution
## x.range:       defined x range (vector with lower and upper)
## plot.line:     bool for drawing the function to the plot
## plot.x.range:  x range for the line to be ploted
## y.lim:         limits for the y axis (vector with lower and upper)
## jump.len:      Jump lenght of the random walk
## norm.const:    normalization constant (applies to the plotted line)
## return.sample: bool for returning the mcmc sample for further use
mcmc <- function(N = 10000, FUN = dnorm, ..., burnin = 0, x0 = 1, discrete = FALSE,
                 x.range = c(-Inf, Inf), plot.line = TRUE, plot.x.range = c(-10,10),
                 y.lim = NULL, jump.len = 1, norm.const = 1, return.sample = FALSE){
                        
  smpl <- c(x0, rep(NA, N-1))
  x = x0
  for (i in 2:N){
    if (discrete){
      if (x > x.range[1] & x < (x.range[2]+jump.len)){
        prop = sample(c(x-jump.len, x + jump.len), prob = c(0.5, 0.5), size = 1)}
      if (x == x.range[2]) {prop = sample(c(x, (x - jump.len)), prob = c(0.5, 0.5), size = 1)}
      if (x == x.range[1]) {prop = sample(c(x, (x + jump.len)), prob = c(0.5, 0.5), size = 1)}
    }
    else{
      while (1) {
        prop = x + runif(1, min = -jump.len, max = jump.len)
        if (prop > x.range[1] & prop < x.range[2]){break()}
      } 
    }
    alpha <- min(1, FUN(prop, ...) / FUN(x, ...))
    x <- sample(c(prop, x), size = 1, replace = FALSE, prob = c(alpha, 1-alpha))
    
    smpl[i] = x
  }
  ret.val = smpl[(burnin+1):length(smpl)]
  typ = "continuous"
  if (discrete){typ = "discrete"}
  lwr = as.character(signif(x.range[1], 3))
  upr = as.character(signif(x.range[2],3 ))
  if (lwr == "-Inf"){lwr = bquote("-\U221E")}
  if (upr == "Inf"){upr = bquote("\U221E")}
  plot.title = paste0(typ, " | x: [", lwr, " ; ", upr, "]")
  if (discrete){
    plot(table(ret.val)/length(ret.val), type = "h", main = plot.title, xlab = "x", ylab = "Density")
    if (plot.line){
      x.sq = plot.x.range[1]:plot.x.range[2]
      lines(x.sq+0.1 , FUN(x.sq, ...)*norm.const, type = "h",col = "red")
    }}
  else{
    hist(ret.val, freq = F, ylim = y.lim, main = plot.title, 
         col = "grey", xlab = "x", breaks = sqrt(length(ret.val)))
    if (plot.line){
      x.sq = seq(plot.x.range[1], plot.x.range[2], by = 0.01)
      lines(x.sq , FUN(x.sq, ...)*norm.const, col = "red")
    }}
  if (return.sample){
    return(ret.val)
  }
}

