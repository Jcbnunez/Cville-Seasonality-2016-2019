### 3. Join analyses for ABC
### 

library(foreach)

root <- "/scratch/yey2sn/Overwintering_ms/20.Connor_sims/Parsed_sims"
files <- system( paste("ls ", root), intern = T )

joint.dat <- foreach(i=1:length(files), .combine = "rbind", .errorhandling = "remove" )%do%{
 
  message(paste(i, "/",length(files)))
  tmp <- get(load( paste( "Parsed_sims/" , files[i], sep = "" ) ))
  
}


