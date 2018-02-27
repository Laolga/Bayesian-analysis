library(mvtnorm)
library(ggmcmc)
library(R2jags)
library(coda)
library(mcmc)
library(xtable)
library(plotly)
setwd("")
#########################################
###############FUNCTIONS#################
#########################################

#initial likelihood for theta and p
L = function(x,y, theta,p) {
  t = x+y
  n = length(x)
  pr = 1
  for (i in 1:n){
    ir = 0
    if (y[i] == 0){
      ir = 1
    }
    pr = pr * ((ir*theta)+((1-theta)*choose(t[i],y[i])*p^(y[i])*(1-p)^(t[i]-y[i])))
  }
  return(pr)
}


#log-likelohood of alpha and delta 
likelihood = function(x,y,alpha,delta){
  return(log(L(x,y,plogis(delta),plogis(alpha))))
}

#log of priors
prior = function(alpha,delta){
  alpha_pr = dnorm(alpha, mean = 0, sd = 10000, log = T)
  delta_pr = dnorm(delta, mean = 0, sd = 10000, log = T)
  return (alpha_pr+delta_pr)
}
#log of posterior
posterior <- function(x,y,alpha,delta){
  return (likelihood(x,y,alpha,delta) + prior(alpha,delta))
}
#METROPOLIS-HASTING algorithm
run_metropolis_MCMC<- function(x,y,startvalue, chain_size, unifa,t,burn_in){
  iterations = chain_size*t+burn_in
  chain = array(dim = c(chain_size,2))
  prev = startvalue
  colnames(chain) = c("alpha","delta")
  chain = as.data.frame(chain)
  ind = 1
  for (i in 1:iterations){
    alpha =  prev[1]
    delta =  prev[2]
    proposal_alpha = alpha + runif(1, min = -unifa, max = unifa)
    proposal_delta = delta + runif(1, min = -unifa, max = unifa)
    
    probab = exp(posterior(x,y,proposal_alpha,proposal_delta) - posterior(x,y,alpha,delta))
    if (runif(1) < probab){
      prev = c(proposal_alpha,proposal_delta)
    }else{
      prev = prev
    }
    if ((i>burn_in)&&(i%%t == 0)){
      chain[ind,]= prev
      ind = ind+1
    }
  }
  return(chain)
}

#Deviance Information Criterion
dic <- function(chain,x,y) {
  n = length(chain[,1])
  lik = rep(NA,n)
  for (i in 1:n){
    lik[i] = likelihood(x,y,chain[i,1],chain[i,2])
  }
  D.bar <- -2*mean(lik)
  theta.bar <- colMeans(chain)
  D.hat <- -2*likelihood(x,y,theta.bar[1],theta.bar[2])
  pD <- D.bar - D.hat
  pV <- var(-2*lik)/2
  list(DIC=pD+D.bar,IC=2*pD+D.bar,pD=pD,pV=pV,Dbar=D.bar,Dhat=D.hat)
}
# Output:
#   list()
#     DIC   : Deviance Information Criterion
#     IC    : Bayesian Predictive Information Criterion
#     pD    : Effective number of parameters (pD = Dbar - Dhat)
#     pV    : Effective number of parameters (pV = var(D)/2)
#     Dbar  : Expected value of the deviance over the posterior
#     Dhat  : Deviance at the mean posterior estimate

#summary of the chain with parameters recovery
chain_sum <-function(chain){
  chain$p = plogis(chain[,1])
  chain$theta = plogis(chain[,3])
  chain$beta = (-1/(chain$p-1))-1
  
  res = list(
    'mean' = colMeans(chain),
    'var' = apply(chain, 2, var),
    'HPD' = HPDinterval(as.mcmc(chain),prob = 0.95),
    't_eff' = effectiveSize(as.mcmc(chain))
  )
  return(res)
}

#latex table for a chain
l_table <- function(chain){
  values = chain_sum(chain) 
  df = data.frame(rbind(t(values$mean),
                        t(values$HPD[,1]),
                        t(values$HPD[,2]),
                        t(values$var),
                        t(values$t_eff)))
  row.names(df) = c('Mean', 'Lower 95%',"Upper 95%",'Var','T_eff')
  table = xtable(df, digits = 4)
  return(table)
}

marginal.likelihood <- function(chain,x,y,log=TRUE) {
  if (class(chain) != "mcmc" & class(chain) != "mcmc.list") {
    stop("chain must be an mcmc or mcmcList object.")
  }
  num.samples = length(chain[,1])
  lik = rep(NA,num.samples)
  for (i in 1:num.samples){
    lik[i] = likelihood(x,y,chain[i,1],chain[i,2])
  }
  
  # get mean parameters
  s <- summary(chain)
  theta.star <- s$statistics[,"Mean"]
  lik.star <- likelihood(x,y,theta.star[1],theta.star[2])
  # get samples from posterior
  g <- sample(1:nrow(chain),num.samples,replace=TRUE)
  theta.g <- chain[g,]
  q.g <- dmvnorm(theta.g,mean=theta.star,sigma=diag(ncol(theta.g)),log=FALSE)
  lik.g = rep(NA,num.samples)
  for (i in 1:num.samples){
    lik.g[i] = likelihood(x,y,theta.g[i,1],theta.g[i,2])
  }
  alpha.g <- sapply(lik.g,function(l) min(1,exp(lik.star-l)))
  # get samples from proposal
  theta.j <- rmvnorm(num.samples,mean=theta.star,sigma=diag(ncol(chain)))
  lik.j = rep(NA,num.samples)
  for (i in 1:num.samples){
    lik.j[i] = likelihood(x,y,theta.j[i,1],theta.j[i,2])
  }
  alpha.j <- sapply(lik.j,function(l) min(1,exp(l-lik.star)))
  pi.hat <- mean(alpha.g*q.g)/mean(alpha.j)
  pi.star <- 0
  
  if (!is.null(prior)) pi.star <- prior(theta.star[1],theta.star[2])
  ln.m <- lik.star + pi.star - log(pi.hat)
  if (! log) ln.m <- exp(ln.m)
  return(ln.m)
}

#bayesian factor 
b_f = function(chain1,chain2){
  bf1 = marginal.likelihood(as.mcmc(chain1),x,y)
  bf2 = marginal.likelihood(as.mcmc(chain2),x,y)
  return(bf1/bf2)
}

batch_var <- function(chain, batch_size = 200){
  n_col = ncol(chain)
  group_split = ceiling(seq_along(chain[,1])/batch_size)
  B = max(group_split)
  b_var = rep(NA, n_col)
  
  for (i in 1:n_col){
    
    samples = split(chain[,i], group_split)
    b_means = lapply(samples, mean)
    i_hat = mean(chain[,i])
    b_var[i] = B * 1 / (batch_size - 1)  * sum((unlist(b_means) - i_hat)^2)
  }
  return (b_var)  
}

#########################################
###############ANALYSYS#################
#########################################
#initial data
x = c(6,9,17,22,7,5,5,14,9,7,9,51)
y = c(5,2,0,0,2,1,0,0,0,0,13,0)

plot(y,col = "blue",type = "b", ylim = c(0,51),lwd = 2, ylab = "PVC", main = "Data representation")
lines(x, col = "orchid",type = "b",lwd = 2)
#lines(t, col = "red",type = "l",lwd = 2)
legend("topleft", c("predrud", "postdrug"), col = c("orchid", "blue"),
       lwd = 2, bty = "n", cex = .9)
#ML estimators
p_seq = t_seq = seq(0.2,0.7,by = 0.003)
MLE = matrix(NA,length(p_seq)*length(t_seq),3)
colnames(MLE) = c("p","theta","MLE")
MLE = data.frame(MLE)
ind = 1
for (i in 1:length(p_seq)){
  for (j in 1:length(t_seq)){
    MLE$MLE[ind] =L(x,y,t_seq[j],p_seq[i])
    MLE$p[ind] =p_seq[i]
    MLE$theta[ind] = t_seq[j]
    ind = ind+1
  }
}

plot_ly(x = MLE$p, y = MLE$theta, z = MLE$MLE,
        marker = list(color = ~mpg, colorscale = c('#FFE1A1', '#683531'), showscale = FALSE)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'p'),
                      yaxis = list(title = 'theta'),
                      zaxis = list(title = 'MLE')),
         annotations = list(
           x = 1.13,
           y = 1.05,
           text = 'Maximum Liklihood',
           xref = 'paper',
           yref = 'paper',
           showarrow = FALSE
         ))
c("theta"=MLE$theta[which.max(MLE$MLE)],"p" = MLE$p[which.max(MLE$MLE)])

#Metropolis_Hasting
chain_1 = run_metropolis_MCMC(x,y,startvalue = c(0,0), chain_size = 30000,unifa =  1.3,t = 1,burn_in = 1000)
с("acceptance rate" = 1-mean(duplicated(chain)))
chain_2 = chain_1[-(1:3000),]
#this one will take long
chain_3 = run_metropolis_MCMC(x,y,startvalue = c(0,0), chain_size = 3800,unifa =  1.3,t = 50,burn_in = 1000)

#jags model
set.seed(123)
dd=list("t" = t,"y" = y, "N"=12)
parameters=c("alpha","delta")
inits=list(alpha = 0, delta = 0)
initial.values=list(inits)
hearts_jags=jags(data=dd,inits=initial.values,parameters.to.save=parameters,model.file="hearts.txt",n.chains=1,n.iter=20000)
print(hearts_jags)
plot(hearts_jags)
traceplot(as.mcmc(chain_3))
loglike_sim=-0.5*hearts_jags$BUGSoutput$sims.matrix[,"deviance"]
chain_j = hearts_jags$BUGSoutput$sims.matrix[,1:2]

#Diagnostics
# to analize chain set varaible "chain" to the name of the mcmc output
#chain =
l_table(chain)
S=ggs(as.mcmc(chain_3))
ggs_compare_partial(S)
ggs_traceplot(S)
ggs_running(S)
ggs_autocorrelation(S,nLags = 100)
ggs_geweke(S)
ggs_caterpillar(S)
heidel.diag(chain)
raftery.diag(chain)
#all together in pdf
#ggmcmc(S)















