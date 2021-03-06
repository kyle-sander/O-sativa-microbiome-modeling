library(deSolve) # https://cran.r-project.org/web/packages/deSolve/vignettes/deSolve.pdf
library(dplyr) # https://www.tidyverse.org/packages 

# this will take the first argument presented in Bash as the # (1-4) of the R instance that it is
args <- commandArgs(trailingOnly = T)
customseed <- as.numeric(args[1])

natcomNumber <- 1000  # how many?
inoculum <- 0.001  # starting abundance
param_mean <- -0.375  # distribution for interaction coefficients
param_sd <- 1
mu_mean <- 0.0002  # distribution for mu (growth rates)
mu_sd <- 0.0002
stabilitythreshold <- 0.1 # determines how strict "timestable" is, used as a percentage of the final abundance
natcomMember <- 5

#differential equation fuction
nativeCommunities <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dX1 <- muX1*X1 + a11*X1*X1 + a12*X1*X2 + a13*X1*X3 + a14*X1*X4 + a15*X1*X5
    dX2 <- muX2*X2 + a22*X2*X2 + a21*X2*X1 + a23*X2*X3 + a24*X2*X4 + a25*X2*X5
    dX3 <- muX3*X3 + a33*X3*X3 + a31*X3*X1 + a32*X3*X2 + a34*X3*X4 + a35*X3*X5
    dX4 <- muX4*X4 + a44*X4*X4 + a41*X4*X1 + a42*X4*X2 + a43*X4*X3 + a45*X4*X5
    dX5 <- muX5*X5 + a55*X5*X5 + a51*X5*X1 + a52*X5*X2 + a53*X5*X3 + a54*X5*X4
    list(c(dX1, dX2, dX3, dX4, dX5))
  })
}

n <- 1
startcomptime <- Sys.time() # for estimating compute time at the end
print(paste0("Begin ", natcomNumber, " Native Community Simulations..."))
set.seed(customseed)

for (n in 1:natcomNumber) {
  print(paste0(n, " of ", natcomNumber))  # shows live progress in console
  
  repeat{
    # starting abundances set by above parameter
    state <- c(X1 = inoculum, X2 = inoculum, X3 = inoculum, X4 = inoculum, X5 = inoculum)  
    
    # interaction parameters and growth rates chosen randomly from distribution supplied above
    natcomParameters <- c(
      muX1 = abs(rnorm(1, mu_mean, sd = mu_sd)), 
      muX2 = abs(rnorm(1, mu_mean, sd = mu_sd)),  
      muX3 = abs(rnorm(1, mu_mean, sd = mu_sd)), 
      muX4 = abs(rnorm(1, mu_mean, sd = mu_sd)), 
      muX5 = abs(rnorm(1, mu_mean, sd = mu_sd)), 
      a11 = abs(rnorm(1, param_mean, sd = param_sd))*-1, 
      a12 = rnorm(1, param_mean, sd = param_sd), 
      a13 = rnorm(1, param_mean, sd = param_sd), 
      a14 = rnorm(1, param_mean, sd = param_sd), 
      a15 = rnorm(1, param_mean, sd = param_sd), 
      a21 = rnorm(1, param_mean, sd = param_sd), 
      a22 = abs(rnorm(1, param_mean, sd = param_sd))*-1, 
      a23 = rnorm(1, param_mean, sd = param_sd), 
      a24 = rnorm(1, param_mean, sd = param_sd), 
      a25 = rnorm(1, param_mean, sd = param_sd), 
      a31 = rnorm(1, param_mean, sd = param_sd), 
      a32 = rnorm(1, param_mean, sd = param_sd), 
      a33 = abs(rnorm(1, param_mean, sd = param_sd))*-1, 
      a34 = rnorm(1, param_mean, sd = param_sd), 
      a35 = rnorm(1, param_mean, sd = param_sd), 
      a41 = rnorm(1, param_mean, sd = param_sd), 
      a42 = rnorm(1, param_mean, sd = param_sd), 
      a43 = rnorm(1, param_mean, sd = param_sd), 
      a44 = abs(rnorm(1, param_mean, sd = param_sd))*-1, 
      a45 = rnorm(1, param_mean, sd = param_sd), 
      a51 = rnorm(1, param_mean, sd = param_sd), 
      a52 = rnorm(1, param_mean, sd = param_sd), 
      a53 = rnorm(1, param_mean, sd = param_sd), 
      a54 = rnorm(1, param_mean, sd = param_sd), 
      a55 = abs(rnorm(1, param_mean, sd = param_sd))*-1)
    
    # 1000 hours in 6 minute intervals
    times <- seq(0, 1000.1, by = 0.1)
    
    # blank matrix for adding abundance data at each timestep as it is generated
    timecourse <- matrix(ncol = 6)
    colnames(timecourse) <- c("time", "X1", "X2", "X3", "X4", "X5")
    
    # new complicated loop for preventing erroneous results:
    z <- 1
    for (z in 1:10001) {
      time1 <- c(times[z], times[z+1]) # only one time step from the 10000
      natcomSimulation <- ode(y = state, times = time1, func = nativeCommunities, parms = natcomParameters)
      # a species whose abuncance drops below 10^-8 is set immediately to 0
      newstate <- as.data.frame(natcomSimulation)[2,] %>% 
        mutate(across(-time, function(x){x <- ifelse(x < 0.00000001, 0, x)}))
      # if the total abuncance rises above 0.1, all species are set to 0, this flags the community for rejection at a later step
      if ((newstate[1,2] + newstate[1,3] + newstate[1,4] + newstate[1,5] + newstate[1,6]) > 0.1) {
        newstate[1,2] <- 0
        newstate[1,3] <- 0
        newstate[1,4] <- 0
        newstate[1,5] <- 0
        newstate[1,6] <- 0
      }
      timecourse <- rbind(timecourse, newstate)
      state <- c(X1 = newstate[1,2], X2 = newstate[1,3], X3 = newstate[1,4], X4 = newstate[1,5], X5 = newstate[1,6])
      z <- z + 1
    }
    
    # calculates rate of change for each species over the final 10 hours
    rocX1 <- abs(timecourse[10000, 2]-timecourse[9900, 2])
    rocX2 <- abs(timecourse[10000, 3]-timecourse[9900, 3])
    rocX3 <- abs(timecourse[10000, 4]-timecourse[9900, 4])
    rocX4 <- abs(timecourse[10000, 5]-timecourse[9900, 5])
    rocX5 <- abs(timecourse[10000, 6]-timecourse[9900, 6])
    
    # rejects a community that has any of the following: a species abundance that is not a number, a total abundance of 0 (indicates a total abundance of 0.1, see above), a species abuncance still changing more than 5% over the last 10 hours
    if (!is.nan(timecourse[10000, 2]) & 
        !is.nan(timecourse[10000, 3]) & 
        !is.nan(timecourse[10000, 4]) & 
        !is.nan(timecourse[10000, 5]) & 
        !is.nan(timecourse[10000, 6]) &
        (timecourse[10000, 2] + timecourse[10000, 3] + timecourse[10000, 4] + timecourse[10000, 5] + timecourse[10000, 6] > 0)) {
      
      # writes the details of a passed community to the output files
      final_states <- matrix(c(paste0("NatCom_", n), timecourse[10000, c(2:6)]), nrow = 1)
      write.table(final_states, file = paste0("UnStable_native_communities_states_", Sys.Date(), ".csv"), sep = ",", col.names = F, row.names = F, append = T)
      
      final_parameters <- matrix(c(paste0("NatCom_", n), natcomParameters), nrow = 1)
      write.table(final_parameters, file = paste0("UnStable_native_communities_parameters_", Sys.Date(), ".csv"), sep = ",", col.names = F, row.names = F, append = T)
      
      final_stability <- matrix(c(paste0("NatCom_", n), rocX1, rocX2, rocX3, rocX4, rocX5), nrow = 1)
      write.table(final_stability, file = paste0("UnStable_native_communities_stability_", Sys.Date(), ".csv"), sep = ",", col.names = F, row.names = F, append = T)
      
      final_timecourse <- timecourse %>% 
        filter(row_number() %% 50 == 0)  # only exports every 5th hour to reduce file size
      final_timecourse <- cbind(rep(paste0("NatCom_", n), nrow(timecourse) / 50), final_timecourse)
      write.table(final_timecourse, file = paste0("UnStable_native_communities_timecourse_", Sys.Date(), ".csv"), sep = ",", col.names = F, row.names = F, append = T)
      
      break
    }
  }
  n <- n + 1
}
# print simulation time
Sys.time() - startcomptime
