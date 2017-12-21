require(reshape2)
require(ggplot2)
require(ggpubr)


# TRACK TOTAL COMPUTING TIME
ptm_total <- proc.time()

# PARAMETERS
total_time <- 1620 # in ms
timesteps <- 400
dt <- total_time/timesteps
time_silence <- 810 # in ms
timestep_silence <- time_silence/total_time*timesteps

cbar_range <- c(0.8,1.4) # min, max cbar values
cap_gamma_range <- c(0.8,2.4)# min and max capital gamma values here
gamma_range <- c(0,4) # min and max gamma

cbar_res <- 20 # how many values do you want to split the cbar range into?
cap_gamma_res <- 20 # how many values do you want to split the cap gamma range into?
gamma_res <- 3 # how many values do you want to split the gamma range into?

num_sims <- 400 # number of simulations for each 'pixel' (parameter combination)

s <- c(-.01,0,.01) # input signals for FIM
M <- 100 # number of individual neurons
tau <- 10 # in ms

xf <- array(0,dim=c(3,num_sims))

# FUNCTIONS
g <- function(x){ # firing rate function
  (tanh(x))
}

L <- function(r){ # LDA conversion function - right now, we're just using the average
  #  -.5*log(v%*%Cpos%*%v) - (((r-mupos)%*%v)^2)/(2*v%*%Cpos%*%v) + .5*log(v%*%Cneg%*%v) + (((r-muneg)%*%v)^2)/(2*v%*%Cneg%*%v)
  mean(r)
}

vt <- function(timestep){ # returns x values at timestep "timestep"
  rowMeans(t(g(x[timestep,,]))) # returns the final averages for each simulation
}

comparis <- function(end_time,init_time){ # evaluate the predictive power, comparing end time to initial time
  results <<- (end_time>0)==(init_time>0)
  return(sum(results)/length(results))
}

stepper <- function(w,range,resolution){ # for iterating across a range of values
  range[1]+(w-1)*((range[2]-range[1])/(resolution-1))
}

nRowMeans<- function(x){ # this function extends rowMeans to arrays with width = 1, instead of just width > 1.
  if(is.null(dim(x))){
    return(x)}
  else {return(rowMeans(x))}
}

#KL divergence using kernel density estimator
density <- function(x, data){
  #density esimator from data using sd=h at new points x
  h <- 1
  n <- length(x) 
  y <- array(0,dim=n)
  for (i in 1:n){
    y[i] <- 1/h * mean(dnorm(x[i], mean=data, sd=h))
  }
  return(y)
}

kl_divergence <- function(data1, data2){
  x = seq(-10,10,length.out = 1000) #initial points when computing kl-divergence
  p <- density(x, data1)
  q <- density(x, data2)
  return(sum(p*log(p/q)))
}

#
# SIMULATIONS
#


mat7a <- array(0,dim=c(cbar_res,cap_gamma_res,gamma_res)) # initialize the results arrays
mat7b <- array(0,dim=c(cbar_res,cap_gamma_res,gamma_res))
mat7c <- array(0,dim=c(cbar_res,cap_gamma_res,gamma_res))
mat7d <- array(0,dim=c(cbar_res,cap_gamma_res,gamma_res))

stepcounter <- 0
avgsteptime <- 0

for (y in 1:gamma_res){
  
  # HETEROGENEITY in connectivity matrix (fixed across cap_gamma, base c simulations)
  gamma <- stepper(y,gamma_range,gamma_res) # fixed gamma, used in building cij (connectivity heterogeneity / noise)
  cijnoise <- matrix(rnorm(M^2,0,gamma),nrow=M,ncol=M)
  cijnoise[lower.tri(cijnoise)] = t(cijnoise)[lower.tri(cijnoise)]   # Make the matrix symmetrical
  diag(cijnoise)<-0
  
  for (w in 1:cbar_res){
    cbar_avg <- stepper(w,cbar_range,cbar_res)
    
    for (z in 1:cap_gamma_res){
      cap_gamma <- stepper(z,cap_gamma_range,cap_gamma_res)
      ptm <- proc.time() # track the time for each iteration
      
      # updated cij matrix:
      cbarbase <- matrix(rep(cbar_avg,M^2),nrow=M,ncol=M)
      diag(cbarbase)<- 0
      cij <- cijnoise + cbarbase
      
      # Initialize for this combination of c, capital gamma:
      lr <- matrix(nrow=num_sims, ncol=M)
      x <- array(0,dim=c(timesteps,M,num_sims))
      xm <- array(0,dim=c(timesteps,M,num_sims))
      xp <- array(0,dim=c(timesteps,M,num_sims))
      
      # Loop through and compute data for each timestep, for (num_sims) trials:
      for(t in 1:timesteps) {
        lastx <- if(t>1) x[(t-1),,] else matrix(0,nrow=M,ncol=num_sims) # set the INITIAL STATE of the neurons here
        dx <- ((ifelse(t<timestep_silence,rep(s[2],M*num_sims),rep(0,M*num_sims)) - lastx + cij%*%g(lastx)/(M-1))*dt + rnorm(M*num_sims,0,sqrt(dt*cap_gamma)))/tau
        x[t,,] <- lastx + dx
      }
      
      for(t in 1:timesteps) { # perturb s in negative direction
        lastx <- if(t>1) xm[(t-1),,] else matrix(0,nrow=M,ncol=num_sims) # set the INITIAL STATE of the neurons here
        dx <- ((ifelse(t<timestep_silence,rep(s[1],M*num_sims),rep(0,M*num_sims)) - lastx + cij%*%g(lastx)/(M-1))*dt + rnorm(M*num_sims,0,sqrt(dt*cap_gamma)))/tau
        xm[t,,] <- lastx + dx
      }
      
      for(t in 1:timesteps) { # perturb s in positive direction
        lastx <- if(t>1) xp[(t-1),,] else matrix(0,nrow=M,ncol=num_sims) # set the INITIAL STATE of the neurons here
        dx <- ((ifelse(t<timestep_silence,rep(s[3],M*num_sims),rep(0,M*num_sims)) - lastx + cij%*%g(lastx)/(M-1))*dt + rnorm(M*num_sims,0,sqrt(dt*cap_gamma)))/tau
        xp[t,,] <- lastx + dx
      }
      
      ############## CALCULATIONS FOR FIGURES ##################
      
      ## 7A ##
      
      mat7a[w,z,y] <- comparis(vt(timesteps),vt(timestep_silence))
      
      ## 7B ##
      
      ranking <- rep(0,M)
      finals <- t(g(x[timesteps,,])) # object contains all the final values for the simulation; each neuron, each sim #!!!!! this is vt(timesteps), I think?
      
      for(zj in 1:M){ # ranks list (descending order) by which neurons are most predictive
        ranking[zj] <- comparis(finals[,zj],vt(timesteps))
      }
      
      finals <- finals[ ,order(ranking, decreasing=T)] # final values ranked by how accurate that particular neuron has been
      
      subsetpred <- rep(0,M)
      for(zi in 1:M){
        subsetpred[zi] <- comparis(nRowMeans(finals[,1:zi]),vt(timestep_silence))
      }
      
      # find the smallest zi for which subsetpred[zi]>=.95*mat7a[w,z,y], which is the predictive power of all the neurons:
      mat7b[w,z,y] <- which(subsetpred>=.95*mat7a[w,z,y])[1]
      
      ## 7C ##
      
      # predictive power as a function of time - figure out the ‘final value’ but for each step, instead of just at the end, and see when it reaches 95% accuracy (across simulations - just exceed 95% of whatever the final accuracy is).
      res <- rep(0,timesteps)
      for(tx in 1:timesteps){
        res[tx] <- comparis(vt(tx),vt(timestep_silence))
      }
      
      # Find the first time at which it goes over .95 final predictive power:
      mat7c[w,z,y] <- which(res>=.95*mat7a[w,z,y])[1]
      
      ## 7D ## FISHER INFORMATION
      
      # average x at final state over all neurons for the three s values:
      xf[1,] <- apply(xm[timesteps,,],c(2),mean) # c(2) means apply across columns
      xf[2,] <- apply(x[timesteps,,],c(2),mean)
      xf[3,] <- apply(xp[timesteps,,],c(2),mean)
      
      FIM <- (kl_divergence(xf[2,], xf[1,]) + kl_divergence(xf[2,], xf[3,]))/(s[2]-s[1])^2
      mat7d[w,z,y] <- FIM # save to corresponding spot in the Figure 7D matrix
  
      
      ##### OUTPUT TO CONSOLE #####
      
      stepcounter <- stepcounter + 1
      
      
      # Console text output:
      cat("\n\nStep ", stepcounter, " of ",cap_gamma_res*gamma_res*cbar_res,"\nCBAR_AVG: ",cbar_avg,"  cap GAMMA: ",cap_gamma,"  gamma: ", gamma)
      #image(mat7a[,,y]) # this updates a preview of the 7a figure as each 'pixel' is computed.
      
      cat("\nFinal predictive ability: ", mat7a[w,z,y])
      cat("\n# neurons for .95 predictive ability: ", mat7b[w,z,y])
      #print(subsetpred)
      cat("\nTimestep of first going over .95 predictive ability: ", mat7c[w,z,y],"\n")
      
      # plot simulation mean trajectories:
      # plot(NULL,xlim=c(min(times),max(times)),ylim=c(-1.25,1.25), main = c("cbar_avg: ",cbar_avg,"cGAM: ",cap_gamma),xlab='time (ms)',ylab='mean neuron state')
      #    title(main = c("cbar_avg: ",cbar_avg,"cGAM: ",cap_gamma), ps=2)
      #for(i in 1:num_sims){points(times,rowMeans(x[,,i]),type="o",col=ifelse(vt(timestep_silence)[i]>0,"red","black"), pch=".")} # prints out results of ea simulation, over time, for each pixel (the first simulation)
      thissimtime <- proc.time()[3] - ptm[3]
      avgsteptime <- (avgsteptime*(stepcounter-1)+thissimtime)/stepcounter
      cat("  Probability correct: ",mat7a[w,z,y],"\n\nTime for this iteration: ", thissimtime,"      [est time remaining: ",(cap_gamma_res*gamma_res*cbar_res-stepcounter)*avgsteptime,"]")
      
    }  # End of loops iterating over cbar, cap gamma values:
  }
}


cat("\nTOTAL runTime: ", proc.time()[3] - ptm_total[3])

##############################
### END OF GENERATING DATA ###
##############################            


# PREPARE DATA FOR PLOTS:
mat7b <- ifelse(mat7a>.6,mat7b,0) # 'black out' pixels in 7b which correspond to final predictive power (7a) below .6

# MAKE PLOTS:
graphs <- function(io){ # Give this function a gamma value, and it will return the plots (Cbar vs capital_gamma)
  f7a <- ggplot(data = melt(mat7a[,,io]), aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() + labs(title="Collective Predictive Power") + scale_x_continuous(name = "C", breaks = c(0,cbar_res/2,cbar_res), labels = c(cbar_range[1],mean(cbar_range),cbar_range[2])) + 
    scale_y_continuous(name = expression(Gamma), breaks = c(0,5,10), labels = c(cap_gamma_range[1],mean(cap_gamma_range),cap_gamma_range[2]))
  
  f7b <- ggplot(data = melt(mat7b[,,io]), aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() + labs(title="# Neurons for 95% Predictive Power") + scale_x_continuous(name = "C", breaks = c(0,cbar_res/2,cbar_res), labels = c(cbar_range[1],mean(cbar_range),cbar_range[2])) + 
    scale_y_continuous(name = expression(Gamma), breaks = c(0,5,10), labels = c(cap_gamma_range[1],mean(cap_gamma_range),cap_gamma_range[2]))
  
  f7c <- ggplot(data = melt(mat7c[,,io]), aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() + labs(title="Timescale of Decision") + scale_x_continuous(name = "C", breaks = c(0,cbar_res/2,cbar_res), labels = c(cbar_range[1],mean(cbar_range),cbar_range[2])) + 
    scale_y_continuous(name = expression(Gamma), breaks = c(0,5,10), labels = c(cap_gamma_range[1],mean(cap_gamma_range),cap_gamma_range[2]))
  
  f7d <- ggplot(data = melt(mat7d[,,io]), aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() + labs(title="FISHER INFORMATION") + scale_x_continuous(name = "C", breaks = c(0,cbar_res/2,cbar_res), labels = c(cbar_range[1],mean(cbar_range),cbar_range[2])) + 
    scale_y_continuous(name = expression(Gamma), breaks = c(0,5,10), labels = c(cap_gamma_range[1],mean(cap_gamma_range),cap_gamma_range[2]))
  
  # paste("Timescale of Decision",stepper(io,gamma_range,gamma_res))
  
  figure <- ggarrange(f7a,f7b,f7c,f7d, labels = c("7a", "7b", "7c","7d"),
                      ncol = 2, nrow = 2, common.legend=FALSE)
  annotate_figure(figure,
                  bottom = text_grob(paste("\u0263: ",round(stepper(io,gamma_range,gamma_res),3)), color = "red", face = "bold", size = 14))
}

# To view graphs, type in graphs(x) where x is the index for your chosen gamma value.