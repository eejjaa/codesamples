#
# PARAMETERS:
#


rep_th <- 0 # minimum c for a parent to reproduce
cost_reproduce <- 0 # cost (in reduction of c) for one reproduction


# INITIALIZE THE WEIGHT OF BIRDS:

ideal_weight_adult <- 80
ideal_weight_chick <- 10

sd_adult <- 8   # std deviation to use in rnorm for adult weight
sd_chick <- 5   # std deviation to use in rnorm for CHICK weight

init_weight <- function(a,idw,sdev){ # give the starting weight(s), sampled from a range:
  return(rnorm(a,mean = idw, sd = sdev))
}

# Weight transformation at maturation:

mw_param <- 1.3 # scaling parameter on the distance from ideal for the mature weight as a function of chick weight.

mat_weight <- function(a,c_weight,sdev){
  ww <- ideal_weight_adult-mw_param*(ideal_weight_chick-c_weight)
  return(rnorm(a,mean = ww, sd = sdev))
}

# OTHER STUFF

age_maturity <- 1 # upon reaching this age, individuals are mature and can reproduce


# # # # # # # # # # #
#     FUNCTIONS     #
# # # # # # # # # # # 

# SURVIVAL functions

# parameters of width for survival function - lower alpha means broader range of weights survive.

alphaC <- .01
alphaP <- .001

Sc <- function(w){
  prob <- 1-alphaC*(ideal_weight_chick-w)^2
  return(ifelse(prob>0,prob,0)) # must be 0 or greater
}

Sp <- function(w){
  prob <- 1-alphaP*(ideal_weight_adult-w)^2
  return(ifelse(prob>0,prob,0)) # must be 0 or greater
}

# DEATH functions

chick_death <- function(w){
  prob <- 1-Sc(w)
  rbinom(1,1,prob)        #      rbinom(1,1,ifelse(prob<1,prob,1)) was used, but since Sc truncates negative values, should be OK
}

parent_death <- function(w){
  prob <- 1-Sp(w)
  rbinom(1,1,prob)
}


# SIGNALING functions

mu <- 0 # metabolic cost to feeding, in units of WEIGHT
V <- 0 # coefficient representing the metabolic cost of begging, in units of WEIGHT
phi <- .2 # expenditure of parent to produce 1 unit of chick weight

xmod_init <- function(a){ # initial values for the measure of exaggeration of need/quality (begging level)
  return(runif(a, min=0, max=.2))
}

ymod_init <- function(a){ # initial values for the adjustment from begging levels to actual resource distribution
  return(runif(a, min=0, max=.2))
}

x <- function(w,xmod){ # solicitation level (x) is some function of offspring state (w)
  return((ideal_weight_chick-w)*(1+xmod)) # alternatives: (1-c)+xmod linearly exaggerate, or (1-c)*(1+xmod), in which exaggeration grows as need grows
}

y <- function(x,ymod,wp){ # resources rendered (in units WEIGHT) is a function only of offspring solicitation level (x)
  y_amt <- x*(1+ymod)
  ifelse(y_amt>=wp,wp,y_amt)
}

# MUTATION functions

mut_sd <- .001 # stdev for a normal distribution of mutation around an xmod or ymod.

mutx <- function(xmod){
#  newval <- xmod + rnorm(1,mean=xmod,sd=mut_sd) # PROBLEMS HERE - getting huge values, etc
#  ifelse(newval>0,newval,0)
#  runif(1,min=xmod-0.01,max=xmod+0.01)
  return(xmod)
}

muty <- function(ymod){
#  newval <- ymod + rnorm(1,mean=ymod,sd=mut_sd)
#  ifelse(newval>0,newval,0)
#  runif(1,min=ymod-0.01,max=ymod+0.01)
  return(ymod)
}

#
# INITIALIZE
#

# SIMULATION PARAMETERS:
init_pop_size <- 900
number_timesteps <- 5000 # number of timesteps (years/generations) to run in the simulation


# Generate initial population:

pop <- data.frame(id=c(1:init_pop_size), type = 1, parent=NA, birthday = 0, deathday = NA, w=init_weight(init_pop_size,ideal_weight_adult,sd_adult), xmod=xmod_init(init_pop_size), ymod=ymod_init(init_pop_size)) # id(or lineage), birthday, f, xmod, ymod

##
## LOOP FOR EACH TIMESTEP
##

identifier <<- init_pop_size + 1
for(t in 1: number_timesteps){
  
  # REPRODUCTION/MATURATION/RESET WEIGHTS:
  
  num_chicks <-0
  
  for(id in pop$id){
    if(pop$type[pop$id==id]==0 && (t-pop$birthday[pop$id==id])>=age_maturity){     # Offspring at or above age_maturity become adults
      pop$type[pop$id==id]<-1
      pop$w[pop$id==id]<-mat_weight(1,pop$w[pop$id==id],sd_adult) # mature weight - should depend on weight at end of 'chick season'
    }
    
    if(pop$w[pop$id==id]>=rep_th && pop$type[pop$id==id]==1){  # Mature organisms reproduce
      pop$w[pop$id==id]<-init_weight(1,ideal_weight_adult,sd_adult)-cost_reproduce # IMPORTANT: no memory of the previous season. If it survived to here, it starts over.
      # also subtract the cost of reproduction from the parent's state
      num_chicks <- num_chicks + 1
      tw <- init_weight(1,ideal_weight_chick,sd_chick) # note that this has nothing to do with parental weight right now
      pop <- rbind(pop,c(identifier,0,id,t,NA,tw,mutx(pop$xmod[pop$id==id]),muty(pop$ymod[pop$id==id]))) # add the offspring row; t= current time, c_init = startup state,
      identifier <<- identifier + 1
    }
  }
  
  
  # SOLICITATION AND FEEDING
  
  for(i in pop$id){
    if(pop$type[pop$id==id]==0 && pop$type[pop$parent[pop$id==id]] == 1){    # if the bird is type 0, AND parent is type 1 (not dead):
      transfer_w <- y(x(pop$w[pop$id==id],pop$xmod[pop$id==id]),pop$ymod[pop$parent[pop$id==id]],pop$w[pop$parent[pop$id==id]]) # make sure strategy calls the PARENT'S ymod!
      pop$w[pop$id==id] <- pop$w[pop$id==id] - V*x(pop$w[pop$id==id],pop$xmod[pop$id==id]) + transfer_w # solicit at cost V*x, with V in WEIGHT units
      tmp <- pop$w[pop$parent[pop$id==id]] - mu - phi*transfer_w # feed at metabolic cost mu (in WEIGHT units), and y*phi
      pop$w[pop$parent[pop$id==id]] <- ifelse(tmp>0,tmp,0)
    }
  }
  
  # GENERIC DENSITY DEPENDENCE: keep the population right around 1000
  popsize <- nrow(pop[pop$type!=2,])
  gdd_temp <- ((popsize - 1000)/popsize)
  gdd <- ifelse(popsize>1000, gdd_temp,0)
  
  for(id in pop$id){
    if(rbinom(1,1,gdd)>0){
      pop$type[pop$id==id]<-2
      pop$deathday[pop$id==id] <- t
    }
  }
  
  
  # DEATHS: Cull based on non-optimal weight.
  
  for(id in pop$id){
    if(pop$type[pop$id==id]==1 && parent_death(pop$w[pop$id==id])==1){  # if parent, check state w/ parent_death(), which may kill;
      pop$type[pop$id==id] <- 2
      pop$deathday[pop$id==id] <- t
    }
    
    if(pop$type[pop$id==id]==0 && chick_death(pop$w[pop$id==id])==1){   # if chick, check state w/ chick_death(), which may kill;
      pop$type[pop$id==id] <- 2
      pop$deathday[pop$id==id] <- t
    }
  }
  
  # REMOVE all of the dead individuals each time, to keep the dataframe smaller:
  pop <- subset(pop, type!=2)

  print(pop)
  cat("Timestep: ", t, "\nPop size: ", nrow(pop))
  print(table(pop$xmod))
  print(table(pop$ymod))
  plot(pop$xmod,pop$ymod) #, xlim=c(0,0.2), ylim=c(0,0.2))
#  readline(prompt="Press [enter] to continue") # pause each step
} # End of each timestep

