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
# POWER ANALYSIS
#



# Sample the data (pull x individuals)

samp <- function(x){  # returns long format of x random animals
    popdex <- sample(7,x,replace =TRUE) # draw from real animals
    #print(popdex)
    scaling <- rgamma(x,3.5,5) #!!!!!!!!!!!!!!!
    #print(scaling)
    
    # turn data into 'long' and add diet:
    df <- melt(vols[popdex,,]*scaling)
    df <- cbind(ifelse(df$Var2<5,0,1),df)
    names(df) <- c("diet", "id","day","band", "volume")
    df$volume <- log(df$volume)
    return(df)
  }

# multiply existing samples by something between, say, .5 and 2, to reduce or exaggerate effect size

corr <- function(d){
  m1 <-  gam(volume ~ s(day) + s(id, bs="re"), data=d) # volumes are already log()
  m0 <-  gam(volume ~    1   + s(id, bs="re"), data=d)
  d$m1 <- m1$fitted.values
  d$vc <- d$volume - d$m1
  pvec <- NA
  for(i in 1:12){
    m1c <-  gam(vc ~ s(day) + s(id, bs="re"), data=d[d$band==i,])
    m0c <-  gam(vc ~    1   + s(id, bs="re"), data=d[d$band==i,])
    pvec[i] <- anova(m1c,m0c,test="Chisq")$`Pr(>Chi)`[2]
  }
  return(pvec) # returns the p-values for each band, after subtracting the overall trend for the dataset.
}

colors <- colorRampPalette(c("yellow", "red","blue"))(20)

w <- function(x){
  plot(NA,main=paste("n = ",x), xlab="Band", ylab="p-val LRT", xlim=c(1,12),ylim=c(0,.06)) # max(a)
  abline(h=.001,col="orange")
  abline(h=.05,col="red")
  for(i in 1:10){
    a <- corr(samp(x))
    points(a,pch=20,col=i)
  }
}

plot(NA, xlab="Band", ylab="p-val LRT", xlim=c(1,12),ylim=c(0,.06)) # max(a)
for(j in 10:20){
  w(j)
}
