require(mgcv)
require(reshape2)
require(ggplot2)
require(ggfortify)
require(stats)


# Function to import from CSV:
get_band <- function(bn){
  data <- as.data.frame(read.csv(paste("band",bn,".csv",sep=""),header=FALSE,stringsAsFactors=FALSE))
  data <- data[9:18,1:44] # trim it down
  data <- (data[3:10,2:44])
  data <- apply(data,FUN=as.double,MARGIN=2)
  a <- c(1:dim(data)[2])
  v <- which(a%%4==1)
  vols[,,bn] <<- data[,v]
  b <- which(a%%4==2)
  bands[,,bn] <<- data[,b]
  w <- which(a%%4==3)
  wts[,,bn] <<- data[,w]
}

# Import:
vols <- array(dim=c(8,11,12)) # 8 individuals, 17 days, 12 protein bands
bands <- array(dim=c(8,11,12)) # 8 individuals, 17 days, 12 protein bands
wts <- array(dim=c(8,11,12)) # 8 individuals, 17 days, 12 protein bands

for(i in 1:12){
  get_band(i)
}


#
# VISUALIZATIONS
#

# plot one animal's band changes over time (color = time):
plot1 <- function(animal){
  #plot.new()
  plot(NULL,xlim=c(1,12),ylim=c(0,max(vols[animal,,],na.rm=T)),xlab="Protein Band",ylab="Expression")
  for(i in 1:11){
    lines(vols[animal,i,], col=colorRampPalette(c("red", "yellow"))(11)[i])
  }
}

# (1) plot all of the animals changing bands over time (time as colors):
plot(NULL,xlim=c(1,12),ylim=c(0,max(vols[,,],na.rm=T)),xlab="Protein Band",ylab="Expression")
for(x in 1:8){
  plot1(x)
}


# (2) heatmaps over time:
hmot <- function(i){
  hm <- ggplot(data = melt(vols[i,,]), aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() + labs(title=paste("Protein Exp over time, Animal",i)) + scale_x_discrete(name = "Day", breaks = c(1:11), 
    labels = c(1:11)) + scale_y_discrete(name = "Protein Band", breaks = c(1:12), labels = c(13,16.5,18.5,19,23,25,34,42,65,75,100,211)) #+
  #scale_fill_manual(values=colors1)
  return(hm)
}

hm_avg <- ggplot(data = melt(apply(vols,c(2,3),mean)), aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + labs(title="Avg Across Animals") + scale_x_discrete(name = "Day", breaks = c(1:11), 
  labels = c(1:11)) + scale_y_discrete(name = "Protein Band", breaks = c(1:12), labels = c(13,16.5,18.5,19,23,25,34,42,65,75,100,211))

hm_var <- ggplot(data = melt(apply(vols,c(2,3),var)), aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + labs(title="Variance") + scale_x_discrete(name = "Day", breaks = c(1:11), 
  labels = c(1:11)) + scale_y_discrete(name = "Protein Band", breaks = c(1:12), labels = c(13,16.5,18.5,19,23,25,34,42,65,75,100,211)) 

ggarrange(hmot(1),hmot(2),hmot(3),hmot(4),hmot(5),hmot(6),hmot(7),hmot(8),hm_avg, labels = c(1:8),
          ncol = 3, nrow = 3, common.legend=TRUE)  

sd <- sqrt(apply(vols,c(2,3),var))

ggplot(data = melt(sd/sqrt(8)), aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + labs(title="STDERR") + scale_x_discrete(name = "Day", breaks = c(1:11), 
  labels = c(1:11)) + scale_y_discrete(name = "Protein Band", breaks = c(1:12), labels = c(13,16.5,18.5,19,23,25,34,42,65,75,100,211)) 



#
# POWER ANALYSIS / MODEL COMPARISON
#

bpv <- function(dfb){
  pvec <- NA
  for(i in 1:12){
    m1<-  gam(log(volume) ~ s(day) + s(id, bs="re"), data=dfb[dfb$band==i,])
    m0<-  gam(log(volume) ~    1   + s(id, bs="re"), data=dfb[dfb$band==i,])
    pvec[i] <- anova(m1,m0,test="Chisq")$`Pr(>Chi)`[2]
  }
  return(pvec)
}

# turn data into 'long' and add diet:
df <- melt(vols)
df <- cbind(ifelse(df$Var2<5,0,1),df)
names(df) <- c("diet", "id","day","band", "volume")


#
# Model comparison approach - randomly sample x individuals from 8 real individuals, re-ID them
# and then build the null vs s(day) models and perform a likelihood ratio test:
#

assess <- function(x){  # returns p-values for each of the 12 bands on a sample of x individual rats
  opts <- c(1,2,3,4,5,6,8)
  popdex <- sample(opts,x,replace =TRUE)
  
  # turn data into 'long' and add diet:
  dft <- melt(vols[popdex,,])
  dft <- cbind(ifelse(dft$Var2<5,0,1),dft)
  names(dft) <- c("diet", "id","day","band", "volume")
  
  return(bpv(dft))
}

pv <- function(sampnum,reps){
  par <- matrix(0,nrow=reps,ncol=12)
  for(i in 1:reps){
    par[i,]<- assess(sampnum)
  }
  
  # visualize:
  plot(colMeans(par),pch=20,col="blue",main="P-val per band", xlab="Band", ylab="p-val LRT", ylim=c(0,.055))
  abline(h=.05,col="red")
  cat("\n# under .05: ", sum(colMeans(par)<.05))
  return(sum(colMeans(par)<.001))
}

# Test across a range of x values and plot the result:

testrange <- c(2:25)
scores <- vector()
for(z in testrange){ # min val = 2
  scores<- c(scores,pv(z,10))
}

plot(x=testrange, y=scores,pch=20,col="orange", main="P-val vs Sample Size", xlab="# rats", ylab="p-val LRT",ylim=c(0,13))
abline(h=12,col="blue")
