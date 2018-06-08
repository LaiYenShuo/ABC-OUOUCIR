rm(list=ls())
### suppose we have data TrueD and suppose we have a normal dist with mean 0, sd = 1
n=10
true.mu <- 10
true.sig <- 5
TrueD <- rnorm(n, mean=true.mu, sd = true.sig)
### we will try approximate ABC MCMC the param of interest are mu and sig^2 suppose prior for mu is normal, prior for sig is uniform
prior <- function(param){
  mu = param[1]
  sig = param[2]
  prior.mu = dnorm(mu, sd=5, log=T)
  prior.sig = dunif(sig, min=0, max=30, log=T)
  return(prior.mu+prior.sig)
}

sumstat <- function(S){
  return(c(mean(S),sd(S)))
}

### we also need proposal function to run abc MCMC
proposalfunction <- function(param){
  return(rnorm(2, mean=param, sd=c((param[1]/2),(param[2]/30))))
}

### to start abc MCMC we need starting value for mu and sigma we run 50000
run_abc_MCMC <- function(startvalue, iterations, TrueD, error=0.05){
  S0 <- sumstat(TrueD)
  chain = array(dim=c(iterations+1,2))
  chain[1,]=startvalue
  for(i in 1:iterations){
    proposal <- proposalfunction(chain[i,])
    ### we generate data using the proposal
    simD <- rnorm(n, mean = proposal[1], sd=proposal[2])
    S1 <- sumstat(simD)
    if(dist(rbind(S0,S1))<error){
      #print(dist(rbind(S0,S1)))
      priorratio <- exp(prior(proposal)-prior(chain[i,]))
      h=min(1,priorratio)
      if(runif(1)<h){
        chain[i+1,]=proposal
        S0 <- S1
        print(c(i,proposal))
      }else{
        chain[i+1,]=chain[i,]
      }
      
    }else{
      chain[i+1,]=chain[i,]
    }
  }
  return(chain)
}
startvalue <- c(10,5)

chain = run_abc_MCMC(startvalue = startvalue, iterations = 100000, TrueD=TrueD, error=0.05)
burnIn=5000
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))
hist(chain[-(1:burnIn),1],nclass=30, main="Posterior of a", xlab="True value = red line")
hist(chain[-(1:burnIn),2],nclass=30, main="Posterior of b", xlab="True value = red line")
