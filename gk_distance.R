### WORK IN PROGRESS ##
# This is a project which compares representations of odors in the glomeruli of the fruit fly brain
# vs the Kenyon cells in the fruit fly brain.


o <- 2 # number of odors
G <- 50 # number of glomeruli
K <- 2000 # number KCs


d <- function(a,b){ # DISTANCE FUNCTION
    #1-sum(a==b)/length(a) #similarity
    1-(a%*%b)/length(a)
}

#initialize these guys
kres <- matrix(data=rep(0,test_range*m_range),nrow=m_range,ncol=test_range)
gres <- matrix(data=rep(0,test_range*m_range),nrow=m_range,ncol=test_range)

m_range <- 10
for(mi in 1:m_range){ # outer loop: glom weights in initial connectivity
  
  ####
  # SET UP THE CONNECTIONS FIRST:
  ####
  
  denscon <- rpois(K,7) # number of connections for each of K KCs
  
  # WEIGHTS of which glomeruli are more likely to be connected to; doing by hand now, trying evenly and over/under rep'd
  glomwts <- sample(10,G,replace=TRUE)
  #glomwts <- abs(rnorm(G,10,mi^2))
  glomwts[1] <- 10*(mi-1)
  cat("\n\nGlom Weights: ",glomwts,"\n")
  # This is the key part to vary across the loop, to try different connectivity patterns
  
  # MAP KCs TO GLOMERULI:
  gk <- matrix(rep(0,G*K),nrow=G,ncol=K)
  
  for (i in 1:K){
    tok <- sample.int(G, size=denscon[i],replace=TRUE,prob=glomwts) # lists which glomeruli KC #i is connected to.
    
    #### NOTE: this was set to replace=FALSE before, maybe that's important!
    
    gk[tok,i]<-1
  }
  
  ####
  # Now test it with different sets of odorants:
  ####
  
  num_samps <- 1000 # number of samples for comparing G distance vs K distance
  test_range <- 100
  
  for (l in 1:test_range){
    condens <- l/test_range # connection density (% glomeruli activated by an odor)
    
    od <- matrix(rbinom(o*G,1,condens),nrow=G,ncol=o) # create an odor-glomerular mapping
    # try different odor densities-of-glomeruli used?
    
    # Map the odors to KC representations:
    m <- matrix(data=0,nrow=K,ncol=o)
    for(i in 1:o){
      active_k <- which(colSums(od[,i]*gk)>1, arr.ind=TRUE) # >1 term reflects 2 or more inputs necessary to activate KC
      m[active_k,i]<-1
    }
    
    gdist <- rep(0,num_samps)
    kdist <- rep(0,num_samps)
    
    for(j in 1:num_samps){
      aj <- sample(o,1)
      bj <- sample(o,1)
      gdist[j] <- d(od[,aj],od[,bj])
      kdist[j] <- d(m[,aj],m[,bj])
    }
    
    cat("\nConnection Density: ",condens, "   G Dist: ",mean(gdist),"   K Dist: ",mean(kdist),"   Ratio K/G: ", mean(kdist)/mean(gdist))
    kres[mi,l]<-mean(kdist)
    gres[mi,l]<-mean(gdist)
    
  }
  
  plot(gres[mi,],title="Red = kenyon, black = Glom", xlab="% of Glomeruli activated per odor",ylab="Avg distance between representations")
  points(kres[mi,],col="red")
  
  # Plot of ratios:
  # In this plot, points below 1 would correspond to more similarity in K reps than in G reps; and vice versa
  # plot(kres/gres, xlab="% of Glomeruli activated per odor", ylab="Ratio of distance between avg K rep and G rep")
}