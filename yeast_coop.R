#require(lattice)


#
# PARAMETERS
#

# BASIC CONFIGURATION:
num_time_steps <- 200 # number of time steps to simulate
env_width <- 30      # width of grid environment
env_height <- 30     # height of grid environment
n <- 100 # Number of yeast organisms to start with
init_glucose <- 5e-7 #6.167497835763118e-7 normally?
init_sucrose <- 2e-7 #6.167497835763118e-7


# CELL LIFE:
#prob_type <- .5 # probability, INITIALLY, that the cell will be type 1 (vs type 0)
divide_energy <- 8e-07 #5e-07 or 06 normally? # threshhold energy needed for an organism to decide to divide
prob_divide <- 1 # probability, given divide_energy, that the organism will divide.
div_range <- 1.0    # max distance away from parent that an offspring can be placed
energy_ratio <- 0.7  # Fraction of the energy that remains with parent organism when splitting
#! cost_of_living <-0 #.1  # Remove this much energy for each timestep
#! old_age_death <- 50   # cells over this many timesteps old will be culled

# ENZYME - SUCROSE/GLUCOSE CONVERSION:

# CONVERSION RATE OF SUCROSE TO GLUCOSE
conv_rate <- function(invJG,S){
  return(ifelse(S>0,(ifelse(((invJG*S)/(10^-8+S))>0,((invJG*S)/(10^-8+S)),0)),0)) # NOTE: previously 10^-4!
}

# INVERTASE ACTIVITY (used for cost, conv_rate evals)
invJG <- function(JG) {
  return(10*(-2.3453*(JG)^2 + 32.609*(JG) + 77.109/(env_width*env_height))) #! Final value is divided over area to reflect specific location. Other coefs not modified!
}
 
# COST TO CONVERT SUCROSE TO GLUCOSE, AS A FUNCTION OF invJG 
cost <- function(invG) {
  #return(0)
  return(10e-07*invG*0.00054*exp(0.000001*invG)) #10e-07 arbitrarily added to reduce cost
}

# EATING
eat_rate <- function(G) {
 return(100*G/(1.52+G)*(1/(env_width*env_height))) # 10*G usually ***-----THESE VALUES ARE FOR 15 MINUTES (1 hr = 40, 15 min = 10), and mmol. 
}

# DIFFUSION
diffcoef<- 3 #!!!!! THIS IS WHERE YOU VARY THE RATE OF DIFFUSION - lower number = more diffusion
diffout <- c(diffcoef,1)
diffin <- c(1,1,1,1,diffcoef) # final digit is the main cell




#
# ORGANISM PROCESSES
#

cell_divide <- function(j,energy_ratio) {
  if(local_orgs[local_orgs$id==j,]$energy > divide_energy && rbinom(1,1,prob_divide)>0){
    local_orgs <<- rbind(local_orgs, local_orgs[local_orgs$id==j,])                     # Create new row for new organism
    popl <<- popl + 1
    local_orgs[(nrow(local_orgs)),]$id <<- popl                                         # New organism ID
    local_orgs[local_orgs$id==popl,]$x <<- local_orgs[local_orgs$id==popl,]$x + runif(1,-div_range,div_range)     # New x-coord
    local_orgs[local_orgs$id==popl,]$y <<- local_orgs[local_orgs$id==popl,]$y + runif(1,-div_range,div_range)     # new y-coord
    local_orgs[local_orgs$id==popl,]$birthday <<- t
    local_orgs[local_orgs$id==popl,]$energy <<- local_orgs[local_orgs$id==j,]$energy*(1-energy_ratio) # New organism's share of energy
    local_orgs[local_orgs$id==j,]$energy <<- local_orgs[local_orgs$id==j,]$energy*energy_ratio # Old organism's share of energy
    }
}
  
enzyme <- function(cellrc,k) {
  # Write to a tracker for organisms$id - how often it applies enzyme, etc (create df for each organism$id, append current timestep for each action)
  if(local_orgs[local_orgs$id==k,]$type==1 && rbinom(1,1,prob_enzyme)>0){
    glucose_layer[cellrc] <<- glucose_layer[cellrc] + (conv_rate(eat_rate(glucose_layer[cellrc]),sucrose_layer[cellrc])) #sucrose_layer[cellrc]*conv_rate is original form
    local_orgs[local_orgs$id==k,]$energy <<- local_orgs[local_orgs$id==k,]$energy - cost(invJG(eat_rate(glucose_layer[cellrc])))
    sucrose_layer[cellrc]<<-sucrose_layer[cellrc] - (conv_rate(eat_rate(glucose_layer[cellrc]),sucrose_layer[cellrc])) #sucrose_layer[cellrc]*conv_rate is original form
  }
}

eat <- function(cellrc,k) {
  # Write to a tracker for organisms$id - how often it eats, etc (create df for each organism$id, append current timestep for each action)
  local_orgs[local_orgs$id==k,]$energy <<- local_orgs[local_orgs$id==k,]$energy + ifelse(glucose_layer[cellrc]>0,eat_rate(glucose_layer[cellrc]),0)
  glucose_layer[cellrc]<<-glucose_layer[cellrc] - ifelse(glucose_layer[cellrc]>0,eat_rate(glucose_layer[cellrc]),0)
}

#
# RESOURCE LAYER FUNCTIONS
#

diffuse <- function(cellrc,layer){
  r <- cellrc[1]
  c <- cellrc[2]
  rN <- ifelse(r==1,env_height,(r-1))
  rS <- ifelse(r==env_height,1,(r+1))
  cW <- ifelse(c==1,env_width,(c-1))
  cE <- ifelse(c==env_width,1,(c+1))
  # base:
  layer[r,c]<-weighted.mean(c(layer[rN,c],layer[rS,c],layer[r,cE],layer[r,cW],layer[r,c]),diffin)
  #cardinals:
  layer[rN,c]<-weighted.mean(c(layer[rN,c],layer[r,c]),diffout)
  layer[rS,c]<-weighted.mean(c(layer[rS,c],layer[r,c]),diffout)
  layer[r,cE]<-weighted.mean(c(layer[r,cE],layer[r,c]),diffout)
  layer[r,cW]<-weighted.mean(c(layer[r,cW],layer[r,c]),diffout)
  
  return(layer)
}

#### START HERE

trackdatalist<-list(NA)

for(typefrac in 0:10){ # start at 0?
    prob_type <- typefrac/10
     
    #################
    # THE MODEL
    #################
    
    
    
    # INITIALIZE ORGANISMS
    organisms <- data.frame(id=1:n, x=runif(n,0,env_width), y=runif(n,0,env_height), type=rbinom(n,1,prob_type), birthday=0, energy=0) # Create a data structure to hold all of the critters
    popl <- n # Initialize the population count
    
    # INITIALIZE RESOURCE LAYERS
    glucose_layer <- matrix(data = init_glucose, nrow=env_height, ncol=env_width) # need to populate these things appropriately
    sucrose_layer <- matrix(data = init_sucrose, nrow=env_height, ncol=env_width)
    index_matrix <- matrix(data=c(1:(env_height*env_width)),nrow=env_height,ncol=env_width,byrow=TRUE) #!!! This needs to be transposed - iterates down cols, not across rows
    
    # INITIALIZE DATA STORAGE
    trackdata <- data.frame(time=0,diff=0,pop0=0,pop1=0,poptot=nrow(organisms),glucose=0,sucrose=0)
    
    #
    grid_list<-sample(nrow(index_matrix)*ncol(index_matrix)) #creates a list of resource cells to iterate through
    
    # THE ACTUAL PROCESS FOR EACH TIME STEP
    for(t in 1:num_time_steps){
      cat("Type Ratio=",prob_type,"\n")
      cat("Time: ",t,"\n")
      cat("Glucose:",mean(glucose_layer)," Sucrose: ",mean(sucrose_layer),"\n")
         
      # LOOP FOR EACH CELL IN GRID
      for(cell in grid_list){
        cellrc <- which(index_matrix==cell,arr.ind=TRUE)
        local_orgs <<- organisms[ceiling(organisms$x)==as.numeric(cellrc[,2]) & ceiling(organisms$y)==as.numeric(cellrc[,1]),]
        # at this point, local_orgs is an array with all the organisms in this cell.
        if(nrow(local_orgs)>0){
          
          # ENZYME
          rand_order <- sample(as.list(local_orgs$id))
          for(org in rand_order){  #!!!
            enzyme(cellrc,org) # apply enzyme for org in cell
          }
   
          # EAT
          rand_order <- sample(as.list(local_orgs$id))
          for(org in rand_order){  #!!!
            eat(cellrc,org) # consume glucose for org in cell
          }
    
          # OTHER ORGANISM FUNCTIONS (Cell division, update strategy, etc)
          for(org in sample(as.list(local_orgs$id))){
            cell_divide(org,energy_ratio) # any organismal reproduction happens here #local_orgs<- if necessary
            update(org) # any bookkeeping functions after each timestep?
          }
          
          # MERGE local_orgs BACK INTO organisms
          for(i in local_orgs$id){
              if(i %in% organisms$id){
                organisms[organisms$id==i,]<-local_orgs[local_orgs$id==i,] # replace rows for existing organisms
              } else {
                organisms<-rbind(organisms,local_orgs[local_orgs$id==i,]) # add rows for new organisms
              }
          }
    #!      cat("CellRC: ",cellrc,"     Number of local_orgs: ", nrow(local_orgs),"\n")
        } # End of IF statement for grid cells which contain organisms
        # Diffuse:
        sucrose_layer <<- diffuse(cellrc,sucrose_layer)
        glucose_layer <<- diffuse(cellrc,glucose_layer)      
    
            
      } # END OF THIS GRID CELL - on to the next
      # BELOW HERE IS STUFF THAT ONLY HAPPENS AFTER EACH TIME STEP
      cat("Type 1 pop: ", nrow(organisms[organisms$type==1,]),"\n")
      cat("Type 0 pop: ", nrow(organisms[organisms$type==0,]),"\n")
      cat("Organisms avg energy:\n")
      print(mean(organisms$energy))
      
      
      # SAVE DATA: trackdata(time    pop0    pop1    glucose sucrose)
    trackdata <- rbind(trackdata, c(time=t,diff=prob_type,pop0=nrow(organisms[organisms$type==0,]),pop1=nrow(organisms[organisms$type==1,]),poptot=nrow(organisms),glucose=sum(glucose_layer),sucrose=sum(sucrose_layer)))


    } # End of for-loop for each time step
    trackdatalist[[length(trackdatalist)+1]]<-trackdata
    #save image:
    png(paste('typeratioplot',prob_type,".png",sep=""), width=1250, height=700)
    plot.new()
    par(mfrow=c(1,2))
    image(sucrose_layer)
    points(organisms[organisms$type==1,]$y/env_width,organisms[organisms$type==1,]$x/env_height,pch=16)
    points(organisms[organisms$type==0,]$y/env_width,organisms[organisms$type==0,]$x/env_height,pch=1)
    image(glucose_layer) #add=TRUE
    points(organisms[organisms$type==1,]$y/env_width,organisms[organisms$type==1,]$x/env_height,pch=16)
    points(organisms[organisms$type==0,]$y/env_width,organisms[organisms$type==0,]$x/env_height,pch=1)
    title(paste("Type Ratio:",prob_type))
    dev.off()
    
} # End of loop running through values of typefrac

##### END #######
  