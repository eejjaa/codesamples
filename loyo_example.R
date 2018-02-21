# Feb 14, 2018

# FUNCTIONS:

pspore <- function(d){ # idea: 95% fall within 1m of individual.
  dnorm(d,mean=0,sd=10) #sd=.5 for 95% within 1 m
}

tree_death <- function(age){
  pd <- pnorm(age,mean=100,sd=35)
  rbinom(1,1,pd)
} # !!! set parameters for deaths here!

mushroom_death <- function(g,v){
  pd <- pnorm(g,mean=30,sd=10)*ifelse(v==-1,2,.5) # !!! cutting prob down if it's not harvested, doubling if it is
  rbinom(1,1,ifelse(pd>1,1,pd))
} # !!! set parameters for deaths here!

makeplots <- function(exp1,exp2){ # makes plot comparing results from experiment 1 and 2
  plot.new()
  par(mfrow=c(2,2))
  plot(exp1$pop,xlab="Population",ylab="",ylim=c(min(exp1$pop,exp2$pop),max(exp1$pop,exp2$pop)))
  points(exp2$pop,col="green")
  plot(exp1$yield,xlab="Yield",ylab="",ylim=c(min(exp1$yield,exp2$yield),max(exp1$yield,exp2$yield)))
  points(exp2$yield,col="green")
  plot(exp1$genvar,xlab="Genetic Diversity",ylab="",pch=20,ylim=c(min(exp1$genvar,exp2$genvar),max(exp1$genvar,exp2$genvar)))
  points(exp2$genvar,col="green",pch=20)
  #plot(exp1$livingages,ylab="Living Ages",xlab="Year")
  #points(exp2$livingages,col="green")
  plot(exp1$deathages,xlab="Death Ages",ylab="",ylim=c(min(exp1$deathages,exp2$deathages),max(exp1$deathages,exp2$deathages)))
  points(exp2$deathages,col="green")
  mtext("Green=Random Harvest, Black=Path Harvest", outer = TRUE) # !!! this isn't working right now
}


#
# PROCESSES:
#

initializ <- function(fr,numtrees){ # fr == rate of trees starting with mushrooms on them
  
  forest <<- data.frame(x=runif(numtrees,0,xdim), y=runif(numtrees,0,ydim), age=abs(floor(rnorm(numtrees,mean=init_tree_mean,sd=init_tree_sd))), 
                        f=rbinom(numtrees,1,fr), g=0, v=0)
  forest$f[forest$f==1] <<- c(1:sum(forest$f)) # numbers (unique ID) for each individual fungus
  forest$g[forest$f>0] <<- abs(floor(rnorm(sum(forest$f>0),mean=30,sd=20))) # !!! make them older - set the age for indivs
  forest$v[forest$f>0 & forest$g>maturity_age] <<- 1 # !!! maturity age? if over 1 year old, then they are viable for reproduction
  genes <<- matrix(0,nrow=length(forest$f[forest$f>0]),ncol=length(forest$f[forest$f>0]))
  diag(genes)<<-1
  deaths <<- data.frame()
}
fragment <- function(fragmentation_type,frag_width,fragtypecontrol){
  
  fw <- .5*frag_width
  #  if(fragmentation_type==0){
  #    # do nothing
  #  }
  
  if(fragmentation_type==4 & fragtypecontrol==TRUE){
    # control quartered:
    forest <<- forest[forest$x>fw & forest$x<xdim-fw,] # trim same amount from edges
    forest <<- forest[forest$y>fw & forest$y<ydim-fw,]
  }
  
  if(fragmentation_type==4 & fragtypecontrol==FALSE){
    # chop into quarters:
    forest <<- forest[forest$x<.5*xdim-fw | forest$x>.5*xdim+fw,] # blocks
    forest <<- forest[forest$y<.5*ydim-fw | forest$y>.5*ydim+fw,]
  }
  
  if(fragmentation_type==16 & fragtypecontrol==TRUE){
    # control sixteenths:
    forest <<- forest[forest$x>3*fw & forest$x<xdim-3*fw,] # trim same amount from edges
    forest <<- forest[forest$y>3*fw & forest$y<ydim-3*fw,]
  }
  
  if(fragmentation_type==16 & fragtypecontrol==FALSE){
    # chop into quarters:
    forest <<- forest[forest$x<.5*xdim-fw | forest$x>.5*xdim+fw,] # blocks
    forest <<- forest[forest$y<.5*ydim-fw | forest$y>.5*ydim+fw,]
    # then into sixteenths:
    forest <<- forest[forest$x<.25*xdim-fw | forest$x>.25*xdim+fw,] # blocks
    forest <<- forest[forest$x<.75*xdim-fw | forest$x>.75*xdim+fw,]
    forest <<- forest[forest$y<.25*ydim-fw | forest$y>.25*ydim+fw,] # blocks
    forest <<- forest[forest$y<.75*ydim-fw | forest$y>.75*ydim+fw,]
  }
} # end of fragment()
propagate <- function(){
  # Compute distance matrix for all trees:
  dm <- as.matrix(dist(forest[,1:2],method="euc",diag=TRUE,upper=TRUE))
  
  # Establish new individuals at viable locations:
  for(k in which(forest$f==0)){ # for each suitable host tree:
    pp <- pspore(dm[k,which(forest$f>0)]) # vec of probs of spores from each possible parent
    probnewindiv <- p.new.indiv # !!! switch to make this dependent on distance
    if(rbinom(1,1,probnewindiv)>0){
      randpars <- sample(forest$f[forest$f>0][pp>0],2,replace=FALSE, prob=pp[pp>0])
      genes<<- rbind(genes,colMeans(rbind(genes[randpars[1],],genes[randpars[2],])))
      
      forest$f[k] <<- max(forest$f)+1
      forest$g[k] <<- 0
      forest$v[k] <<- 0 # starts out not being viable (needs to be over 5 yr old in this version)
    }
  }
}
update_fungi <- function(){
  
  # AGING:
  forest$g[forest$f>0]<<-forest$g[forest$f>0]+1
  
  # NATURAL DEATHS:
  m <- mapply(mushroom_death, g=forest$g, v=forest$v) # probabilities of death based on age, harvest status
  tokill <- which(m==1) # which mushrooms are being killed off (tree IDs)
  if(length(tokill)>0){
    ages <- cbind(forest$f[tokill],forest$g[tokill])
    deaths <<- rbind(deaths,ages)
  }
  forest$v[tokill]<<-0 # clear out the dead ones
  forest$f[tokill]<<-0
  forest$g[tokill]<<-0
  
  # MATURATION: (return to viability):
  forest$v[forest$g>5] <<- 1 # as they age, they become viable
  forest$v[forest$v==-1] <<- 1 # harvested previous year, bouncing back.
  
}
update_trees <- function(){
  # Kill old trees (and the fungi on them):
  z <- sapply(forest$age,tree_death) # probabilities of death based on age in function tree_death()
  dyf <- forest[which(z==1),] # trees marked for death
  
  if(sum(dyf$f>0)>0){ # if there are any mushrooms on the dying trees,
    ages <- cbind(dyf$f[dyf$f>0],dyf$g[dyf$f>0]) # assemble their IDs and ages
    deaths <<- rbind(deaths,ages)
  }
  forest <<- forest[which(z==0),] # then remove the dead trees from the forest
  num_dead_trees <- sum(z)
  
  # Add new trees:
  nt <- rpois(1,num_dead_trees) # fancier old version: abs(round(rnorm(1,mean=tree_reproduce_frac*nrow(forest),sd=6))) # number of new trees to add
  #cat("\n\n\n nt: ",nt)
  forest <<- rbind(forest, data.frame(x=runif(nt,0,1000), y=runif(nt,0,1000), age=rep(0,nt), 
                                      f=rep(0,nt), g=rep(0,nt),v=rep(0,nt)))
  # VARIATION - FOR SHIFTING FOREST (rightward over time:
  # forest <<- rbind(forest, data.frame(x=runif(nt,i*10,1000), y=runif(nt,i*10,1000), age=rep(0,nt), 
  #                                     f=rep(0,nt), g=rep(0,nt),v=rep(0,nt)))
  
  # Aging of trees:
  forest$age <<- forest$age+1
}
update_trees_cc <- function(){
  # Kill old trees (and the fungi on them):
  z <- sapply(forest$age,tree_death) # probabilities of death based on age in function tree_death()
  z <- c(z,sapply(forest$x,cc_death)) # uses cc_death to kill off stuff that is farther left.
  dyf <- forest[which(z==1),] # trees marked for death
  
  if(sum(dyf$f>0)>0){ # if there are any mushrooms on the dying trees,
    ages <- cbind(dyf$f[dyf$f>0],dyf$g[dyf$f>0]) # assemble their IDs and ages
    deaths <<- rbind(deaths,ages)
  }
  forest <<- forest[which(z==0),] #then remove the dead trees from the forest
  num_dead_trees <- sum(z)
  
  # Add new trees:
  nt <- rpois(1,num_dead_trees) # fancier old version: abs(round(rnorm(1,mean=tree_reproduce_frac*nrow(forest),sd=6))) # number of new trees to add
  #cat("\n\n\n nt: ",nt)
  forest <<- rbind(forest, data.frame(x=runif(nt,0,1000), y=runif(nt,0,1000), age=rep(0,nt), 
                                      f=rep(0,nt), g=rep(0,nt),v=rep(0,nt)))
  # VARIATION - FOR SHIFTING FOREST (rightward over time:
  # forest <<- rbind(forest, data.frame(x=runif(nt,i*10,1000), y=runif(nt,i*10,1000), age=rep(0,nt), 
  #                                     f=rep(0,nt), g=rep(0,nt),v=rep(0,nt)))
  
  # Aging of trees:
  forest$age <<- forest$age+1
}
harvest_random <- function(harvestnumber){
  viable <<- sum(forest$v==1)
  s <- min(rpois(1,harvestnumber),sum(forest$v==1))
  m <- which(forest$v==1) # which trees contain mature (harvestable) loyo
  # h <- which(rbinom(length(m),1,s)==1) # take a fraction of the population (old way)
  h <- sample(m,s) #choose around s of those and collect (with poisson lambda = s)
  forest$v[h] <<- -1 # turn those individuals' viability off
  return(length(h)) # this function also returns the yield, the number of individuals harvested.
}
harvest_path <- function(harvestnumber){ # this slows down over iterations! why?
  viable <<- sum(forest$v==1)
  fruits <- forest[forest$v==1,]
  tdf <- as.matrix(dist(fruits[,1:2],method="euc",diag=TRUE,upper=TRUE))
  
  s <- min(rpois(1,harvestnumber),sum(forest$v==1))
  
  # choose a mushroom to start with:
  x <- sample(nrow(fruits),1)
  closdist <- 0
  fruits$v[x]<--1
  
  for(j in 2:s){ # choose approximately s (poisson lambda = s)
    # harvest x;
    #ds<-cbind(tdf[,x],fruits$v)
    # find the closest one:
    closdist <- min(tdf[,x][fruits$v==1])
    # assign it to be x:
    x <- which(as.vector(tdf[,x])==closdist,arr.ind=TRUE)
    fruits$v[x]<--1
    #print(dists)
    #print(x)
    #print(fruits[x,1:2])
    #print(fruits$v)
    points(fruits[fruits$v==-1,],pch=19,cex=.5,col="chocolate3")
    points(forest[forest$v==-1,],pch=19,cex=.5,col="chocolate3")
  }
  forest$v[forest$v==1] <<- fruits$v
  return(sum(fruits$v==-1)) # this function also returns the yield, the number of individuals harvested.
}
track <- function(h,i){
  
  living_ages <- forest$g[forest$g>0]
  livingagesvec[h,i] <<- mean(living_ages)
  deathagesvec[h,i] <<- mean(deaths[,2])
  genvar[h,i] <<- var(colMeans(genes[-deaths[,1],]))
  
  popvec[h,i] <<- sum(forest$f>0)
  yieldvec[h,i] <<- yield
  
  # Console Output:
  cat("\n\n ",i)#,"\nPop: ",pop[i],"  mean tree age: ",mean(forest$age),
  cat("\ntrees: ",nrow(forest),"  mushrooms: ",sum(forest$f>0)," viable: ",viable," YIELD: ",yield)
  #cat("\nProcTime: ",proc.time()[3]-ptm)
}
update_plot <- function(sf){
  plot(forest$x,forest$y,xlim=c(0,xdim),ylim=c(0,ydim),pch=8,col=ifelse(forest$age>mean(forest$age),"forestgreen","dark green"), 
       cex=ifelse(forest$age>mean(forest$age),2*sf,2.3*sf), lwd=2*sf, axes=FALSE)
  
  points(forest[forest$f>0,],pch=20,cex=1.8*sf,col="goldenrod")
  points(forest[forest$v==1,],pch=1,cex=1.5*sf,col="coral3",lwd=sf)
  points(forest[forest$v==1,],pch=20,cex=.7*sf,col="coral3")
  points(forest[forest$v==-1,],pch=1,cex=2*sf,col="yellow",lwd=sf) #coral1
}

#
# SIMULATIONS:
#

# Initial conditions:
p.new.indiv <- .01
init_tree_mean <- 60
init_tree_sd <- 5

tree_reproduce_frac <- .05
maturity_age <- 5

# Constant across simulations:
xdim <- 1000
ydim <- 1000

wrapper <- function(plot_token,number_simulations,numyears, numtrees,     tree_reproduce_frac,init_prevalence_fungus,  harvest_strategy,harvest_number, fragmentation_type,frag_width,fragtypecontrol){
  
  # create dataframes to hold number_simulations instances of pop, yieldvec:
  popvec <<- matrix(NA,nrow=number_simulations,ncol=numyears)
  livingagesvec <<- matrix(NA,nrow=number_simulations,ncol=numyears)
  deathagesvec <<- matrix(NA,nrow=number_simulations,ncol=numyears)
  yieldvec <<- matrix(NA,nrow=number_simulations,ncol=numyears)
  genvar <<- matrix(NA,nrow=number_simulations,ncol=numyears)
  
  for(h in 1:number_simulations){ # Start of each SIM
    cat("\nSim #",h," of ",number_simulations)
    ptm <- proc.time()[3]
    initializ(init_prevalence_fungus,numtrees)
    fragment(fragmentation_type,frag_width,fragtypecontrol)
    if(plot_token == TRUE){ update_plot(1) }
    
    for(i in 1:numyears){ # Start of each YEAR
      i<<-i
      if(harvest_strategy==1){
        yield <<- harvest_random(harvest_number)
      }
      if(harvest_strategy==2){
        yield <<- harvest_path(harvest_number)
      }
      if(plot_token == TRUE){ update_plot(1) }
      propagate()
      update_trees()
      update_fungi()
      track(h,i)
      #readline(prompt="Press [enter] to continue") # DEBUG
    } # end of each 'year'
    cat("\nTotal time for Sim: ",proc.time()[3]-ptm)
  } # end of each SIM
  # Across all simulations:
  output <-list(pop=colMeans(popvec), yield=colMeans(yieldvec), livingages=colMeans(livingagesvec), deathages=colMeans(deathagesvec), genvar=colMeans(genvar)) # each of these dataframes/vectors is an element in the output list.
  return(output)
} # end of wrapper()


# Example calls: comparing two harvesting strategies, with all else constant

pf4ctrl <- wrapper(plot_token=FALSE,number_simulations=10,numyears=1000, numtrees=1000,     tree_reproduce_frac=.05,init_prevalence_fungus=.2,
                   harvest_strategy=2,harvest_number=100, fragmentation_type=4,frag_width=10,fragtypecontrol=TRUE)

rf4ctrl <- wrapper(plot_token=FALSE,number_simulations=10,numyears=1000, numtrees=1000,     tree_reproduce_frac=.05,init_prevalence_fungus=.2,
                   harvest_strategy=1,harvest_number=100, fragmentation_type=4,frag_width=10,fragtypecontrol=TRUE)