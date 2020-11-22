library(rootSolve)
library(dplyr)

wrkdir <- "~/SynComEvo/TimeStable_natcoms_201120/"

# this will take the first argument presented in Bash as the number of the PSC to be evolved
args <- commandArgs(trailingOnly = T)
PSC_number <- as.numeric(args[1]) #file number name where KillerComm and Single Killer files are stored

# leave this alone
num_of_comms <- 1000
num_of_species <- 11
# assigns initial conditions, in this case all species start at same abundance
start_amt <- 0.001
num_of_Common_KillerCommParams <- 80 #number of parameters to cycle through and, separately, set to zero and assess performance

# gLV function for single killer simulations
NonLinearGLV6<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    dX1 <- muX1*X1 + a11*X1*X1 + a12*X1*X2 + a13*X1*X3 + a14*X1*X4 + a15*X1*X5 + a16*X1*X6
    dX2 <- muX2*X2 + a22*X2*X2 + a21*X2*X1 + a23*X2*X3 + a24*X2*X4 + a25*X2*X5 + a26*X2*X6
    dX3 <- muX3*X3 + a33*X3*X3 + a31*X3*X1 + a32*X3*X2 + a34*X3*X4 + a35*X3*X5 + a36*X3*X6
    dX4 <- muX4*X4 + a44*X4*X4 + a41*X4*X1 + a42*X4*X2 + a43*X4*X3 + a45*X4*X5 + a46*X4*X6
    dX5 <- muX5*X5 + a55*X5*X5 + a51*X5*X1 + a52*X5*X2 + a53*X5*X3 + a54*X5*X4 + a56*X5*X6
    dX6 <- muX6*X6 + a66*X6*X6 + a61*X6*X1 + a62*X6*X2 + a63*X6*X3 + a64*X6*X4 + a65*X6*X5
    list(c(dX1,dX2,dX3,dX4,dX5,dX6))
  })
}

#gLV function for probiotic support community simulations
NonLinearGLV11 <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dX1 <- muX1*X1 + a11*X1*X1 + a12*X1*X2 + a13*X1*X3 + a14*X1*X4 + a15*X1*X5 + a16*X1*X6 + a17*X1*X7 + a18*X1*X8 + a19*X1*X9 + a110*X1*X10 + a1_11*X1*X11
    dX2 <- muX2*X2 + a21*X2*X1 + a22*X2*X2 + a23*X2*X3 + a24*X2*X4 + a25*X2*X5 + a26*X2*X6 + a27*X2*X7 + a28*X2*X8 + a29*X2*X9 + a210*X2*X10 + a211*X2*X11
    dX3 <- muX3*X3 + a31*X3*X1 + a32*X3*X2 + a33*X3*X3 + a34*X3*X4 + a35*X3*X5 + a36*X3*X6 + a37*X3*X7 + a38*X3*X8 + a39*X3*X9 + a310*X3*X10 + a311*X3*X11
    dX4 <- muX4*X4 + a41*X4*X1 + a42*X4*X2 + a43*X4*X3 + a44*X4*X4 + a45*X4*X5 + a46*X4*X6 + a47*X4*X7 + a48*X4*X8 + a49*X4*X9 + a410*X4*X10 + a411*X4*X11
    dX5 <- muX5*X5 + a51*X5*X1 + a52*X5*X2 + a53*X5*X3 + a54*X5*X4 + a55*X5*X5 + a56*X5*X6 + a57*X5*X7 + a58*X5*X8 + a59*X5*X9 + a510*X5*X10 + a511*X5*X11
    dX6 <- muX6*X6 + a61*X6*X1 + a62*X6*X2 + a63*X6*X3 + a64*X6*X4 + a65*X6*X5 + a66*X6*X6 + a67*X6*X7 + a68*X6*X8 + a69*X6*X9 + a610*X6*X10 + a611*X6*X11
    dX7 <- muX7*X7 + a71*X7*X1 + a72*X7*X2 + a73*X7*X3 + a74*X7*X4 + a75*X7*X5 + a76*X7*X6 + a77*X7*X7 + a78*X7*X8 + a79*X7*X9 + a710*X7*X10 + a711*X7*X11
    dX8 <- muX8*X8 + a81*X8*X1 + a82*X8*X2 + a83*X8*X3 + a84*X8*X4 + a85*X8*X5 + a86*X8*X6 + a87*X8*X7 + a88*X8*X8 + a89*X8*X9 + a810*X8*X10 + a811*X8*X11
    dX9 <- muX9*X9 + a91*X9*X1 + a92*X9*X2 + a93*X9*X3 + a94*X9*X4 + a95*X9*X5 + a96*X9*X6 + a97*X9*X7 + a98*X9*X8 + a99*X9*X9 + a910*X9*X10 + a911*X9*X11
    dX10 <- muX10*X10 + a101*X10*X1 + a102*X10*X2 + a103*X10*X3 + a104*X10*X4 + a105*X10*X5 + a106*X10*X6 + a107*X10*X7 + a108*X10*X8 + a109*X10*X9 + a1010*X10*X10 + a1011*X10*X11
    dX11 <- muX11*X11 + a11_1*X11*X1 + a112*X11*X2 + a113*X11*X3 + a114*X11*X4 + a115*X11*X5 + a116*X11*X6 + a117*X11*X7 + a118*X11*X8 + a119*X11*X9 + a1110*X11*X10 + a1111*X11*X11
    list(c(dX1, dX2, dX3, dX4, dX5, dX6, dX7, dX8, dX9, dX10, dX11))
  })
}

#reading in parameters from native community generation
csvFilePath <- paste0(wrkdir, list.files(wrkdir)[grep("native_communities_parameters_", list.files(wrkdir))])
NativeCommunities_stode_final_parameters <- read.csv(file = csvFilePath, header = FALSE)

#reading in states from native community generation
csvFilePath1 <- paste0(wrkdir, list.files(wrkdir)[grep("native_communities_states_", list.files(wrkdir))])
NativeCommunities_stode_final_states <- read.csv(file = csvFilePath1, header = FALSE)

# initialize matrix for tracking important parms
importantParms <- matrix(nrow = 0, ncol = 11)
colnames(importantParms) <- c("Parameter Number being omitted/set to 0", "Original Parameter Value", "# of NC's target total abundance > 0.001", "# of NC target species Diff > 0 (total abundance)", "#num of NC with all negative eigenvalues", "% performance when parameter omitted", "Parameter", "# of NC target species [0-single target inoc] -> [>0 - PSC inoc]", "% Performance Difference [0-single target inoc] -> [>0 - PSC inoc]", "# of NC target species [>0 - single target inoc] -> [>>0 - PSC inoc] ", "% Performance Difference [>0 - single target inoc] -> [>>0 - PSC inoc]")

# start loop that is repeated for each round of evolution
z <- 1
repeat{ 
  # sets variables used for comparison in the first round to avoid confusion
  if (z == 1) {
    topPSC <- PSC_number
  }
  oldcheck <- importantParms[, 1]
  
  print("parameter sweeping...")
  
  #create directory for new data
  dir.create(paste0(wrkdir, "evolution_of_PSC_", PSC_number, "/Round_", z, "_data"))
  
  #reading in selected PSC parameters 
  csvFilePath2 <- paste0(wrkdir, "evolution_of_PSC_", PSC_number, "/Round_", (z-1), "_data", "/KillerComm_", topPSC, "_stode_final_parameters.csv")
  KillerComm_stode_highPerformingParams <- read.csv(file = csvFilePath2, header = T)
  
  #reading in selected PSC states
  csvFilePath3 <- paste0(wrkdir, "evolution_of_PSC_", PSC_number, "/Round_", (z-1), "_data", "/KillerComm_", topPSC, "_stode_final_states.csv")
  KillerComm_stode_highPerformingStates <- read.csv(file = csvFilePath3, header = FALSE)
  
  #reading in selected PSC eigs
  csvFilePath4 <- paste0(wrkdir, "evolution_of_PSC_", PSC_number, "/Round_", (z-1), "_data", "/KillerComm_", topPSC, "_stode_final_eigs.csv")
  KillerComm_stode_highPerformingEigs <- read.csv(file = csvFilePath4, header = FALSE)
  
  #reading in single target species inoculation final states
  csvFilePath5 <- paste0(wrkdir, "evolution_of_PSC_", PSC_number, "/Round_", (z-1), "_data", "/SingleKiller_stode_final_states.csv")
  SingleKillerFinStates <- read.csv(file = csvFilePath5, header = FALSE)
  
  #reading in single target species inoculation final parameters
  csvFilePath6 <- paste0(wrkdir, "evolution_of_PSC_", PSC_number, "/Round_", (z-1), "_data", "/SingleKiller_stode_final_parameters.csv")
  SingleKiller_stode_final_parameters <- read.csv(file = csvFilePath6, header = T)
  
  #reading in metrics file from the previous round
  csvFilePath7 <- paste0(wrkdir, "evolution_of_PSC_", PSC_number, "/Round_", (z-1), "_data", "/", list.files(paste0(wrkdir, "evolution_of_PSC_", PSC_number, "/Round_", (z-1), "_data"))[grep("KillerCommMetrics_file_number_", list.files(paste0(wrkdir, "evolution_of_PSC_", PSC_number, "/Round_", (z-1), "_data")))])
  KillerCommMetrics <- read.csv(file = csvFilePath7, stringsAsFactors = F)
  
  # initialize matrix for storing parameter sweep data
  ParameterSensitivity <- matrix(nrow = num_of_Common_KillerCommParams + 1, ncol = 11)
  col_names <- c("Parameter Number being omitted/set to 0", "Original Parameter Value", "# of NC's target total abundance > 0.001", "# of NC target species Diff > 0 (total abundance)", "#num of NC with all negative eigenvalues", "% performance when parameter omitted", "Parameter", "# of NC target species [0-single target inoc] -> [>0 - PSC inoc]", "% Performance Difference [0-single target inoc] -> [>0 - PSC inoc]", "# of NC target species [>0 - single target inoc] -> [>>0 - PSC inoc] ", "% Performance Difference [>0 - single target inoc] -> [>>0 - PSC inoc]")
  colnames(ParameterSensitivity) <- col_names
  
  #Import KillerComm_stode_final_states as well to calculate baseline, no parameter change and insert into first row of table
  ParameterSensitivity[1, 1] <- paste("baseline - no parameters omitted")  #This is identifier for scenario when no parameters have been augmented/omitted, the 'baseline' - 
 
   #use this matrix position as row label to indicate as such - instead of "0"
  ParameterSensitivity[1, 3] <- sum(KillerComm_stode_highPerformingStates[, 6] > 0.001, na.rm = TRUE) #this is calculating which communities are >0.001 right now
  
  #want to shift this to be diffs > 0 as needed parameter - its cool this is next column, need column headers 
  KillerCommTotalDiffs <- KillerComm_stode_highPerformingStates[, 6] - SingleKillerFinStates[, 6]
  ParameterSensitivity[1, 4] <- sum(KillerCommTotalDiffs > 0, na.rm = TRUE)  # of Native Comms KillerComm target species Diff > 0 (total)

  num_of_neg_eigs <- rowSums(KillerComm_stode_highPerformingEigs < 0)  
  ParameterSensitivity[1, 5] <- length(which(num_of_neg_eigs == num_of_species))  # num of (-) eigenvalues for single target species inoculation
  
  ParameterSensitivity[1, 8] <- as.numeric(KillerCommMetrics[topPSC, 21]) #reading in [0]->[>0] for baseline, from KillerCommMetrics
  
  ParameterSensitivity[1, 10] <- as.numeric(KillerCommMetrics[topPSC, 22]) #reading in [>0]->[>>0] for baseline, from KillerCommMetrics
  
  #begin KillerComm community generation, reset loop counters
  n <- 1  # loop counter for number of native + killer combos
  p <- 1  # loop counter for number of KillerComm Common Parameters individually set to zero and reassessed
  
  # begins loop that will be repeated, setting each individual parameter of CommonKillerCommParams, separately, 
  # to zero and then re-assessing support community performance against all 1000 native communities
  for (p in 1:num_of_Common_KillerCommParams) {
    
    # initialize matrices for storing single parameter sweep data
    zero_parameters_KillerComm <- matrix(ncol = (num_of_species * num_of_species + num_of_species), nrow = num_of_comms)
    zero_states_KillerComm <- matrix(ncol = num_of_species, nrow = num_of_comms)
    zero_eigs_KillerComm <- matrix(ncol = num_of_species, nrow = num_of_comms)
    
    #import parameters from KillerComms that are common to all Native Comms
    CommonKillerCommParams <- matrix(nrow = 1, ncol = num_of_Common_KillerCommParams)
    CommonKillerCommParams[1,1] <- KillerComm_stode_highPerformingParams$a17[[1000]]
    CommonKillerCommParams[1,2] <- KillerComm_stode_highPerformingParams$a18[[1000]]
    CommonKillerCommParams[1,3] <- KillerComm_stode_highPerformingParams$a19[[1000]]
    CommonKillerCommParams[1,4] <- KillerComm_stode_highPerformingParams$a110[[1000]]
    CommonKillerCommParams[1,5] <- KillerComm_stode_highPerformingParams$a1_11[[1000]]
    CommonKillerCommParams[1,6] <- KillerComm_stode_highPerformingParams$a27[[1000]]
    CommonKillerCommParams[1,7] <- KillerComm_stode_highPerformingParams$a28[[1000]]
    CommonKillerCommParams[1,8] <- KillerComm_stode_highPerformingParams$a29[[1000]]
    CommonKillerCommParams[1,9] <- KillerComm_stode_highPerformingParams$a210[[1000]]
    CommonKillerCommParams[1,10] <- KillerComm_stode_highPerformingParams$a211[[1000]]
    CommonKillerCommParams[1,11] <- KillerComm_stode_highPerformingParams$a37[[1000]]
    CommonKillerCommParams[1,12] <- KillerComm_stode_highPerformingParams$a38[[1000]]
    CommonKillerCommParams[1,13] <- KillerComm_stode_highPerformingParams$a39[[1000]]
    CommonKillerCommParams[1,14] <- KillerComm_stode_highPerformingParams$a310[[1000]]
    CommonKillerCommParams[1,15] <- KillerComm_stode_highPerformingParams$a311[[1000]]
    CommonKillerCommParams[1,16] <- KillerComm_stode_highPerformingParams$a47[[1000]]
    CommonKillerCommParams[1,17] <- KillerComm_stode_highPerformingParams$a48[[1000]]
    CommonKillerCommParams[1,18] <- KillerComm_stode_highPerformingParams$a49[[1000]]
    CommonKillerCommParams[1,19] <- KillerComm_stode_highPerformingParams$a410[[1000]]
    CommonKillerCommParams[1,20] <- KillerComm_stode_highPerformingParams$a411[[1000]]
    CommonKillerCommParams[1,21] <- KillerComm_stode_highPerformingParams$a57[[1000]]
    CommonKillerCommParams[1,22] <- KillerComm_stode_highPerformingParams$a58[[1000]]
    CommonKillerCommParams[1,23] <- KillerComm_stode_highPerformingParams$a59[[1000]]
    CommonKillerCommParams[1,24] <- KillerComm_stode_highPerformingParams$a510[[1000]]
    CommonKillerCommParams[1,25] <- KillerComm_stode_highPerformingParams$a511[[1000]]
    CommonKillerCommParams[1,26] <- KillerComm_stode_highPerformingParams$a67[[1000]]
    CommonKillerCommParams[1,27] <- KillerComm_stode_highPerformingParams$a68[[1000]]
    CommonKillerCommParams[1,28] <- KillerComm_stode_highPerformingParams$a69[[1000]]
    CommonKillerCommParams[1,29] <- KillerComm_stode_highPerformingParams$a610[[1000]]
    CommonKillerCommParams[1,30] <- KillerComm_stode_highPerformingParams$a611[[1000]]
    CommonKillerCommParams[1,31] <- KillerComm_stode_highPerformingParams$a71[[1000]]
    CommonKillerCommParams[1,32] <- KillerComm_stode_highPerformingParams$a72[[1000]]
    CommonKillerCommParams[1,33] <- KillerComm_stode_highPerformingParams$a73[[1000]]
    CommonKillerCommParams[1,34] <- KillerComm_stode_highPerformingParams$a74[[1000]]
    CommonKillerCommParams[1,35] <- KillerComm_stode_highPerformingParams$a75[[1000]]
    CommonKillerCommParams[1,36] <- KillerComm_stode_highPerformingParams$a76[[1000]]
    CommonKillerCommParams[1,37] <- KillerComm_stode_highPerformingParams$a78[[1000]]
    CommonKillerCommParams[1,38] <- KillerComm_stode_highPerformingParams$a79[[1000]]
    CommonKillerCommParams[1,39] <- KillerComm_stode_highPerformingParams$a710[[1000]]
    CommonKillerCommParams[1,40] <- KillerComm_stode_highPerformingParams$a711[[1000]]
    CommonKillerCommParams[1,41] <- KillerComm_stode_highPerformingParams$a81[[1000]]
    CommonKillerCommParams[1,42] <- KillerComm_stode_highPerformingParams$a82[[1000]]
    CommonKillerCommParams[1,43] <- KillerComm_stode_highPerformingParams$a83[[1000]]
    CommonKillerCommParams[1,44] <- KillerComm_stode_highPerformingParams$a84[[1000]]
    CommonKillerCommParams[1,45] <- KillerComm_stode_highPerformingParams$a85[[1000]]
    CommonKillerCommParams[1,46] <- KillerComm_stode_highPerformingParams$a86[[1000]]
    CommonKillerCommParams[1,47] <- KillerComm_stode_highPerformingParams$a87[[1000]]
    CommonKillerCommParams[1,48] <- KillerComm_stode_highPerformingParams$a89[[1000]]
    CommonKillerCommParams[1,49] <- KillerComm_stode_highPerformingParams$a810[[1000]]
    CommonKillerCommParams[1,50] <- KillerComm_stode_highPerformingParams$a811[[1000]]
    CommonKillerCommParams[1,51] <- KillerComm_stode_highPerformingParams$a91[[1000]]
    CommonKillerCommParams[1,52] <- KillerComm_stode_highPerformingParams$a92[[1000]]
    CommonKillerCommParams[1,53] <- KillerComm_stode_highPerformingParams$a93[[1000]]
    CommonKillerCommParams[1,54] <- KillerComm_stode_highPerformingParams$a94[[1000]]
    CommonKillerCommParams[1,55] <- KillerComm_stode_highPerformingParams$a95[[1000]]
    CommonKillerCommParams[1,56] <- KillerComm_stode_highPerformingParams$a96[[1000]]
    CommonKillerCommParams[1,57] <- KillerComm_stode_highPerformingParams$a97[[1000]]
    CommonKillerCommParams[1,58] <- KillerComm_stode_highPerformingParams$a98[[1000]]
    CommonKillerCommParams[1,59] <- KillerComm_stode_highPerformingParams$a910[[1000]]
    CommonKillerCommParams[1,60] <- KillerComm_stode_highPerformingParams$a911[[1000]]
    CommonKillerCommParams[1,61] <- KillerComm_stode_highPerformingParams$a101[[1000]]
    CommonKillerCommParams[1,62] <- KillerComm_stode_highPerformingParams$a102[[1000]]
    CommonKillerCommParams[1,63] <- KillerComm_stode_highPerformingParams$a103[[1000]]
    CommonKillerCommParams[1,64] <- KillerComm_stode_highPerformingParams$a104[[1000]]
    CommonKillerCommParams[1,65] <- KillerComm_stode_highPerformingParams$a105[[1000]]
    CommonKillerCommParams[1,66] <- KillerComm_stode_highPerformingParams$a106[[1000]]
    CommonKillerCommParams[1,67] <- KillerComm_stode_highPerformingParams$a107[[1000]]
    CommonKillerCommParams[1,68] <- KillerComm_stode_highPerformingParams$a108[[1000]]
    CommonKillerCommParams[1,69] <- KillerComm_stode_highPerformingParams$a109[[1000]]
    CommonKillerCommParams[1,70] <- KillerComm_stode_highPerformingParams$a1011[[1000]]
    CommonKillerCommParams[1,71] <- KillerComm_stode_highPerformingParams$a11_1[[1000]]
    CommonKillerCommParams[1,72] <- KillerComm_stode_highPerformingParams$a112[[1000]]
    CommonKillerCommParams[1,73] <- KillerComm_stode_highPerformingParams$a113[[1000]]
    CommonKillerCommParams[1,74] <- KillerComm_stode_highPerformingParams$a114[[1000]]
    CommonKillerCommParams[1,75] <- KillerComm_stode_highPerformingParams$a115[[1000]]
    CommonKillerCommParams[1,76] <- KillerComm_stode_highPerformingParams$a116[[1000]]
    CommonKillerCommParams[1,77] <- KillerComm_stode_highPerformingParams$a117[[1000]]
    CommonKillerCommParams[1,78] <- KillerComm_stode_highPerformingParams$a118[[1000]]
    CommonKillerCommParams[1,79] <- KillerComm_stode_highPerformingParams$a119[[1000]]
    CommonKillerCommParams[1,80] <- KillerComm_stode_highPerformingParams$a1110[[1000]]
  
    #set parameter number p from CommonKillerCommParams equal to zero
    ParameterSensitivity[p + 1, 1] <- p  #denotes which index of p is currently running
    ParameterSensitivity[p + 1, 2] <- CommonKillerCommParams[1, p]  #denotes original 
    CommonKillerCommParams[1, p] <- 0  

    #combine Common KillerComm Params and other KillerComm params in single vector
    onezeroKCparms <- c(
      muX7 = KillerComm_stode_highPerformingParams$muX7[[1000]],
      muX8 = KillerComm_stode_highPerformingParams$muX8[[1000]],
      muX9 = KillerComm_stode_highPerformingParams$muX9[[1000]],
      muX10 = KillerComm_stode_highPerformingParams$muX10[[1000]],
      muX11 = KillerComm_stode_highPerformingParams$muX11[[1000]],
      a17 = CommonKillerCommParams[1,1],
      a18 = CommonKillerCommParams[1,2],
      a19 = CommonKillerCommParams[1,3],
      a110 = CommonKillerCommParams[1,4],
      a1_11 = CommonKillerCommParams[1,5],
      a27 = CommonKillerCommParams[1,6],
      a28 = CommonKillerCommParams[1,7],
      a29 = CommonKillerCommParams[1,8],
      a210 = CommonKillerCommParams[1,9],
      a211 = CommonKillerCommParams[1,10],
      a37 = CommonKillerCommParams[1,11],
      a38 = CommonKillerCommParams[1,12],
      a39 = CommonKillerCommParams[1,13],
      a310 = CommonKillerCommParams[1,14],
      a311 = CommonKillerCommParams[1,15],
      a47 = CommonKillerCommParams[1,16],
      a48 = CommonKillerCommParams[1,17],
      a49 = CommonKillerCommParams[1,18],
      a410 = CommonKillerCommParams[1,19],
      a411 = CommonKillerCommParams[1,20],
      a57 = CommonKillerCommParams[1,21],
      a58 = CommonKillerCommParams[1,22],
      a59 = CommonKillerCommParams[1,23],
      a510 = CommonKillerCommParams[1,24],
      a511 = CommonKillerCommParams[1,25],
      a67 = CommonKillerCommParams[1,26],
      a68 = CommonKillerCommParams[1,27],
      a69 = CommonKillerCommParams[1,28],
      a610 = CommonKillerCommParams[1,29],
      a611 = CommonKillerCommParams[1,30],
      a71 = CommonKillerCommParams[1,31],
      a72 = CommonKillerCommParams[1,32],
      a73 = CommonKillerCommParams[1,33],
      a74 = CommonKillerCommParams[1,34],
      a75 = CommonKillerCommParams[1,35],
      a76 = CommonKillerCommParams[1,36],
      a77 = KillerComm_stode_highPerformingParams$a77[[1000]],
      a78 = CommonKillerCommParams[1,37],
      a79 = CommonKillerCommParams[1,38],
      a710 = CommonKillerCommParams[1,39],
      a711 = CommonKillerCommParams[1,40],
      a81 = CommonKillerCommParams[1,41],
      a82 = CommonKillerCommParams[1,42],
      a83 = CommonKillerCommParams[1,43],
      a84 = CommonKillerCommParams[1,44],
      a85 = CommonKillerCommParams[1,45],
      a86 = CommonKillerCommParams[1,46],
      a87 = CommonKillerCommParams[1,47],
      a88 = KillerComm_stode_highPerformingParams$a88[[1000]],
      a89 = CommonKillerCommParams[1,48],
      a810 = CommonKillerCommParams[1,49],
      a811 = CommonKillerCommParams[1,50],
      a91 = CommonKillerCommParams[1,51],
      a92 = CommonKillerCommParams[1,52],
      a93 = CommonKillerCommParams[1,53],
      a94 = CommonKillerCommParams[1,54],
      a95 = CommonKillerCommParams[1,55],
      a96 = CommonKillerCommParams[1,56],
      a97 = CommonKillerCommParams[1,57],
      a98 = CommonKillerCommParams[1,58],
      a99 = KillerComm_stode_highPerformingParams$a99[[1000]],
      a910 = CommonKillerCommParams[1,59],
      a911 = CommonKillerCommParams[1,60],
      a101 = CommonKillerCommParams[1,61],
      a102 = CommonKillerCommParams[1,62],
      a103 = CommonKillerCommParams[1,63],
      a104 = CommonKillerCommParams[1,64],
      a105 = CommonKillerCommParams[1,65],
      a106 = CommonKillerCommParams[1,66],
      a107 = CommonKillerCommParams[1,67],
      a108 = CommonKillerCommParams[1,68],
      a109 = CommonKillerCommParams[1,69],
      a1010 = KillerComm_stode_highPerformingParams$a1010[[1000]],
      a1011 = CommonKillerCommParams[1,70],
      a11_1 = CommonKillerCommParams[1,71],
      a112 = CommonKillerCommParams[1,72],
      a113 = CommonKillerCommParams[1,73],
      a114 = CommonKillerCommParams[1,74],
      a115 = CommonKillerCommParams[1,75],
      a116 = CommonKillerCommParams[1,76],
      a117 = CommonKillerCommParams[1,77],
      a118 = CommonKillerCommParams[1,78],
      a119 = CommonKillerCommParams[1,79],
      a1110 = CommonKillerCommParams[1,80],
      a1111 = KillerComm_stode_highPerformingParams$a1111[[1000]])
    
    # begins loop that will be repeated for each combo of killercomm + native comm
    for (n in 1:num_of_comms) {
      # reading in parameters from single killer model and native communities
      NCandSKparms <- c(
        muX1 = NativeCommunities_stode_final_parameters[n, 1],
        muX2 = NativeCommunities_stode_final_parameters[n, 2],
        muX3 = NativeCommunities_stode_final_parameters[n, 3],
        muX4 = NativeCommunities_stode_final_parameters[n, 4],
        muX5 = NativeCommunities_stode_final_parameters[n, 5],
        a11 = NativeCommunities_stode_final_parameters[n, 6],
        a12 = NativeCommunities_stode_final_parameters[n, 7],
        a13 = NativeCommunities_stode_final_parameters[n, 8],
        a14 = NativeCommunities_stode_final_parameters[n, 9],
        a15 = NativeCommunities_stode_final_parameters[n, 10],
        a21 = NativeCommunities_stode_final_parameters[n, 11],
        a22 = NativeCommunities_stode_final_parameters[n, 12],
        a23 = NativeCommunities_stode_final_parameters[n, 13],
        a24 = NativeCommunities_stode_final_parameters[n, 14],
        a25 = NativeCommunities_stode_final_parameters[n, 15],
        a31 = NativeCommunities_stode_final_parameters[n, 16],
        a32 = NativeCommunities_stode_final_parameters[n, 17],
        a33 = NativeCommunities_stode_final_parameters[n, 18],
        a34 = NativeCommunities_stode_final_parameters[n, 19],
        a35 = NativeCommunities_stode_final_parameters[n, 20],
        a41 = NativeCommunities_stode_final_parameters[n, 21],
        a42 = NativeCommunities_stode_final_parameters[n, 22],
        a43 = NativeCommunities_stode_final_parameters[n, 23],
        a44 = NativeCommunities_stode_final_parameters[n, 24],
        a45 = NativeCommunities_stode_final_parameters[n, 25],
        a51 = NativeCommunities_stode_final_parameters[n, 26],
        a52 = NativeCommunities_stode_final_parameters[n, 27],
        a53 = NativeCommunities_stode_final_parameters[n, 28],
        a54 = NativeCommunities_stode_final_parameters[n, 29],
        a55 = NativeCommunities_stode_final_parameters[n, 30],
        muX6 = SingleKiller_stode_final_parameters[n, 31],
        a16 = SingleKiller_stode_final_parameters[n, 32],
        a61 = SingleKiller_stode_final_parameters[n, 33],
        a26 = SingleKiller_stode_final_parameters[n, 34],
        a62 = SingleKiller_stode_final_parameters[n, 35],
        a36 = SingleKiller_stode_final_parameters[n, 36],
        a63 = SingleKiller_stode_final_parameters[n, 37],
        a46 = SingleKiller_stode_final_parameters[n, 38],
        a64 = SingleKiller_stode_final_parameters[n, 39],
        a56 = SingleKiller_stode_final_parameters[n, 40],
        a65 = SingleKiller_stode_final_parameters[n, 41],
        a66 = SingleKiller_stode_final_parameters[n, 42])
      
      # sets starting abundance to the final states obtained from the previous models and adds killercomm abundances
      state <- c(X1 = NativeCommunities_stode_final_states[n, 1],
                 X2 = NativeCommunities_stode_final_states[n, 2],
                 X3 = NativeCommunities_stode_final_states[n, 3],
                 X4 = NativeCommunities_stode_final_states[n, 4],
                 X5 = NativeCommunities_stode_final_states[n, 5],
                 X6 = start_amt,
                 X7 = start_amt,
                 X8 = start_amt,
                 X9 = start_amt,
                 X10 = start_amt,
                 X11 = start_amt)

      # combine native + single killer parameters with the killer comm parameters into a single vector for the non-linear solver
      parametersKillerComm <- c(NCandSKparms, onezeroKCparms)
      
      # the non-linear solver
      out_4 <- stode(y = state, fun = NonLinearGLV11, parms = parametersKillerComm, pos = TRUE)
      
      # calculate jacobian for eigenvalue evaluation of stability
      jac <- jacobian.full(y = c(out_4$y[1], out_4$y[2], out_4$y[3], out_4$y[4], out_4$y[5], out_4$y[6], out_4$y[7], out_4$y[8], out_4$y[9], out_4$y[10], out_4$y[11]), func = NonLinearGLV11, parms = parametersKillerComm)
      eig <- eigen(jac)
      
      # adds final states, params, and eigs to matrices. Only the states matrix is ever used. The others may help troubleshooting.
      zero_states_KillerComm[n, ] = out_4$y
      zero_parameters_KillerComm[n, ] <- parametersKillerComm
      zero_eigs_KillerComm[n, ] <- Re(eig$values)
      
      n = n + 1
    }
    
    # num of native comms target total abundance > 0.001
    ParameterSensitivity[p + 1, 3] <- sum(zero_states_KillerComm[, 6] > 0.001, na.rm = TRUE) 
    
    # num of native comms target total abundance diff from single killer > 0
    KillerCommTotalDiffs <- zero_states_KillerComm[, 6] - SingleKillerFinStates[, 6]
    ParameterSensitivity[p + 1, 4] <- sum(KillerCommTotalDiffs > 0, na.rm = TRUE)# of Native Comms KillerComm target species Diff > 0 (total)
    
    # num of native comms with all Re(eigs) < 0  
    num_of_neg_eigs <- rowSums(zero_eigs_KillerComm < 0)  
    ParameterSensitivity[p + 1, 5] <- length(which(num_of_neg_eigs == num_of_species)) #num of NC with all negative eigenvalues
    
    # % of performance difference (NC Target Species>0 when param omitted/NC Target Species>0 when no params omitted)*100
    #NOT WORKING ??
    ParameterSensitivity[p + 1, 6] <- as.numeric(ParameterSensitivity[p + 1, 4]) / as.numeric(ParameterSensitivity[1, 4]) * 100
    
    # num of NC [0 in single species inoc] -> [>0 in probiotic support community inoc]
      #1. [Single Killer Fin States - species 6] AND [final States Killer Comm - species 6] into same matrix, two columns
      perform_metrics_table <- matrix(nrow = num_of_comms, ncol = 5)
      perform_metrics_table[, 1] <- zero_states_KillerComm[, 6]
      perform_metrics_table[, 2] <- SingleKillerFinStates[, 6]
      #2. column 3, give a 1 if single killer Fin state = 0
      perform_metrics_table[, 3] <- ifelse(perform_metrics_table[, 2] == 0, 1, 0)
      #3. column 4, give a 1 if [killercomm-singleinoc]>0
      perform_metrics_table[, 4] <- ifelse(perform_metrics_table[, 1] > perform_metrics_table[, 2], 1, 0)
      #4. column 5, add two columns
      perform_metrics_table[, 5] <- perform_metrics_table[, 3] + perform_metrics_table[, 4]
      #5. NC of [0] -> [>0] is sums(column 5 = 2) - number of rows whose column #5 is = 2
      NCs_zero_to_pos <- sum(perform_metrics_table[, 5] == 2)
      #6. Put this # of NC inot ParameterSensitivity[p+1,7]
    ParameterSensitivity[p + 1, 8] <- NCs_zero_to_pos
    
    # calculating % difference in [0 in single species inoc] -> [>0 in probiotic support community inoc] with parameter omitted
    pct_diff_zeroToPos <- (NCs_zero_to_pos / as.numeric(ParameterSensitivity[1, 8])) * 100 
    ParameterSensitivity[p + 1, 9] <- pct_diff_zeroToPos
    
    # of NC [>0 in single species inoc] -> [>>0 in probiotic support community inoc] 
      #1. [Single Killer Fin States - species 6] AND [final States Killer Comm - species 6] into same matrix, two columns
      perform_metrics_table2 <- matrix(nrow = num_of_comms, ncol = 5)
      perform_metrics_table2[, 1] <- zero_states_KillerComm[, 6] #final state of target species from PSC
      perform_metrics_table2[, 2] <- SingleKillerFinStates[, 6] #final states of target species from single target species inoculation
      #2. column 3, give a 1 if single killer Fin state > 0
      perform_metrics_table2[, 3] <- ifelse(perform_metrics_table2[, 2] > 0, 1, 0)
      #3. column 4, give a 1 if [killercomm-singleinoc]>0
      perform_metrics_table2[, 4] <- ifelse(perform_metrics_table2[, 1] > perform_metrics_table2[, 2], 1, 0)
      #4. column 5, add two columns
      perform_metrics_table2[, 5] <- perform_metrics_table2[, 3] + perform_metrics_table2[, 4]
      #5. NC of [0] -> [>0] is sums(column 5 = 2) - number of rows whose column #5 is = 2
      NCs_pos_to_greater <- sum(perform_metrics_table2[, 5] == 2)
      #6. Put this # of NC inot ParameterSensitivity[p+1,8]
    ParameterSensitivity[p + 1, 10] <- NCs_pos_to_greater
   
    # calculating % difference in [>0 in single species inoc] -> [>>0 in probiotic support community inoc] with parameter omitted
    pct_diff_posToGreater <- (NCs_pos_to_greater / as.numeric(ParameterSensitivity[1, 10])) * 100 
    ParameterSensitivity[p + 1, 11] <- pct_diff_posToGreater
    
    rm(perform_metrics_table)
    rm(perform_metrics_table2)

      if (p %% 10 == 0){
      print(p)
    }

    p = p + 1
    n = 1
  }
  
  # fill in parameter names manually to parametersensitivity table
  ParameterSensitivity[2, 7] <- paste("a17")
  ParameterSensitivity[3, 7] <- paste("a18")
  ParameterSensitivity[4, 7] <- paste("a19")
  ParameterSensitivity[5, 7] <- paste("a110")
  ParameterSensitivity[6, 7] <- paste("a1_11")
  
  ParameterSensitivity[7, 7] <- paste("a27")
  ParameterSensitivity[8, 7] <- paste("a28")
  ParameterSensitivity[9, 7] <- paste("a29")
  ParameterSensitivity[10, 7] <- paste("a210")
  ParameterSensitivity[11, 7] <- paste("a211")
  
  ParameterSensitivity[12, 7] <- paste("a37")
  ParameterSensitivity[13, 7] <- paste("a38")
  ParameterSensitivity[14, 7] <- paste("a39")
  ParameterSensitivity[15, 7] <- paste("a310")
  ParameterSensitivity[16, 7] <- paste("a311")
  
  ParameterSensitivity[17, 7] <- paste("a47")
  ParameterSensitivity[18, 7] <- paste("a48")
  ParameterSensitivity[19, 7] <- paste("a49")
  ParameterSensitivity[20, 7] <- paste("a410")
  ParameterSensitivity[21, 7] <- paste("a411")
  
  ParameterSensitivity[22, 7] <- paste("a57")
  ParameterSensitivity[23, 7] <- paste("a58")
  ParameterSensitivity[24, 7] <- paste("a59")
  ParameterSensitivity[25, 7] <- paste("a510")
  ParameterSensitivity[26, 7] <- paste("a511")
  
  ParameterSensitivity[27, 7] <- paste("a67")
  ParameterSensitivity[28, 7] <- paste("a68")
  ParameterSensitivity[29, 7] <- paste("a69")
  ParameterSensitivity[30, 7] <- paste("a610")
  ParameterSensitivity[31, 7] <- paste("a611")
  
  ParameterSensitivity[32, 7] <- paste("a71")
  ParameterSensitivity[33, 7] <- paste("a72")
  ParameterSensitivity[34, 7] <- paste("a73")
  ParameterSensitivity[35, 7] <- paste("a74")
  ParameterSensitivity[36, 7] <- paste("a75")
  ParameterSensitivity[37, 7] <- paste("a76")
  ParameterSensitivity[38, 7] <- paste("a78")
  ParameterSensitivity[39, 7] <- paste("a79")
  ParameterSensitivity[40, 7] <- paste("a710")
  ParameterSensitivity[41, 7] <- paste("a711")
  
  ParameterSensitivity[42, 7] <- paste("a81")
  ParameterSensitivity[43, 7] <- paste("a82")
  ParameterSensitivity[44, 7] <- paste("a83")
  ParameterSensitivity[45, 7] <- paste("a84")
  ParameterSensitivity[46, 7] <- paste("a85")
  ParameterSensitivity[47, 7] <- paste("a86")
  ParameterSensitivity[48, 7] <- paste("a87")
  ParameterSensitivity[49, 7] <- paste("a89")
  ParameterSensitivity[50, 7] <- paste("a810")
  ParameterSensitivity[51, 7] <- paste("a811")
  
  ParameterSensitivity[52, 7] <- paste("a91")
  ParameterSensitivity[53, 7] <- paste("a92")
  ParameterSensitivity[54, 7] <- paste("a93")
  ParameterSensitivity[55, 7] <- paste("a94")
  ParameterSensitivity[56, 7] <- paste("a95")
  ParameterSensitivity[57, 7] <- paste("a96")
  ParameterSensitivity[58, 7] <- paste("a97")
  ParameterSensitivity[59, 7] <- paste("a98")
  ParameterSensitivity[60, 7] <- paste("a910")
  ParameterSensitivity[61, 7] <- paste("a911")
  
  ParameterSensitivity[62, 7] <- paste("a101")
  ParameterSensitivity[63, 7] <- paste("a102")
  ParameterSensitivity[64, 7] <- paste("a103")
  ParameterSensitivity[65, 7] <- paste("a104")
  ParameterSensitivity[66, 7] <- paste("a105")
  ParameterSensitivity[67, 7] <- paste("a106")
  ParameterSensitivity[68, 7] <- paste("a107")
  ParameterSensitivity[69, 7] <- paste("a108")
  ParameterSensitivity[70, 7] <- paste("a109")
  ParameterSensitivity[71, 7] <- paste("a1011")
  
  ParameterSensitivity[72, 7] <- paste("a11_1")
  ParameterSensitivity[73, 7] <- paste("a112")
  ParameterSensitivity[74, 7] <- paste("a113")
  ParameterSensitivity[75, 7] <- paste("a114")
  ParameterSensitivity[76, 7] <- paste("a115")
  ParameterSensitivity[77, 7] <- paste("a116")
  ParameterSensitivity[78, 7] <- paste("a117")
  ParameterSensitivity[79, 7] <- paste("a118")
  ParameterSensitivity[80, 7] <- paste("a119")
  ParameterSensitivity[81, 7] <- paste("a1110")
  
  # write full Parameter Sensitivity table into file
  writeFilePath <- paste0(wrkdir, "evolution_of_PSC_", PSC_number, "/Round_", z, "_data", "/ParameterSensitivity_", z, "_PSC_", topPSC, ".csv")
  write.table(ParameterSensitivity, file = writeFilePath, sep = ",", row.names = FALSE)
  
  # trim Parameter Sensitivity to keep only those with an impact greater than -10%
  newParms <- as.data.frame(ParameterSensitivity, stringsAsFactors = F) %>% 
    mutate(`% Performance Difference [>0 - single target inoc] -> [>>0 - PSC inoc]` = as.numeric(`% Performance Difference [>0 - single target inoc] -> [>>0 - PSC inoc]`)) %>% 
    filter(`% Performance Difference [>0 - single target inoc] -> [>>0 - PSC inoc]` < 90)

  
  importantParms <- rbind(importantParms, newParms)
  
  parmcheck <- newParms[, 1]
  
  # write partial Parameter Sensitivity table into file
  writeFilePath <- paste0(wrkdir, "evolution_of_PSC_", PSC_number, "/selected_ParameterSensitivity_log.csv")
  write.table(importantParms, file = writeFilePath, sep = ",", row.names = FALSE)
  
  # used for comparing progress to previous rounds
  parmcount <- nrow(newParms)
  
  #initializing variables for new killer comm simulations
  num_of_KillerComm <- 2000 #set back to 2000
  # mean and standard deviation for interaction coefficients  
  param_mean <- (-0.375)
  param_sd <- 1
  # mean and standard deviation for growth rates
  mu_mean <- 0.0002
  mu_sd <- 0.0002
  
  final_parameters_singlekiller <- matrix(ncol = 42, nrow = num_of_comms)
  colnames(final_parameters_singlekiller) <- c("muX1",
                                               "muX2",
                                               "muX3",
                                               "muX4",
                                               "muX5",
                                               "a11",
                                               "a12",
                                               "a13",
                                               "a14",
                                               "a15",
                                               "a21",
                                               "a22",
                                               "a23",
                                               "a24",
                                               "a25",
                                               "a31",
                                               "a32",
                                               "a33",
                                               "a34",
                                               "a35",
                                               "a41",
                                               "a42",
                                               "a43",
                                               "a44",
                                               "a45",
                                               "a51",
                                               "a52",
                                               "a53",
                                               "a54",
                                               "a55",
                                               "muX6",
                                               "a16",
                                               "a61",
                                               "a26",
                                               "a62",
                                               "a36",
                                               "a63",
                                               "a46",
                                               "a64",
                                               "a56",
                                               "a65",
                                               "a66")
 
  final_states_singlekiller <- matrix(ncol = 6, nrow = num_of_comms)
  
  final_eigs_singlekiller <- matrix(ncol=6,nrow=num_of_comms)
  n <- 1 #begin loop counter 

  SKparamsForSKreassess <- SingleKiller_stode_final_parameters
  
  #begin single killer community generation
  for (n in 1:num_of_comms) {
    # import native community and single killer parameters to vector
    parametersSingleKiller <- c(
      muX1 = NativeCommunities_stode_final_parameters[n, 1],
      muX2 = NativeCommunities_stode_final_parameters[n, 2],
      muX3 = NativeCommunities_stode_final_parameters[n, 3],
      muX4 = NativeCommunities_stode_final_parameters[n, 4],
      muX5 = NativeCommunities_stode_final_parameters[n, 5],
      a11 = NativeCommunities_stode_final_parameters[n, 6],
      a12 = NativeCommunities_stode_final_parameters[n, 7],
      a13 = NativeCommunities_stode_final_parameters[n, 8],
      a14 = NativeCommunities_stode_final_parameters[n, 9],
      a15 = NativeCommunities_stode_final_parameters[n, 10],
      a21 = NativeCommunities_stode_final_parameters[n, 11],
      a22 = NativeCommunities_stode_final_parameters[n, 12],
      a23 = NativeCommunities_stode_final_parameters[n, 13],
      a24 = NativeCommunities_stode_final_parameters[n, 14],
      a25 = NativeCommunities_stode_final_parameters[n, 15],
      a31 = NativeCommunities_stode_final_parameters[n, 16],
      a32 = NativeCommunities_stode_final_parameters[n, 17],
      a33 = NativeCommunities_stode_final_parameters[n, 18],
      a34 = NativeCommunities_stode_final_parameters[n, 19],
      a35 = NativeCommunities_stode_final_parameters[n, 20],
      a41 = NativeCommunities_stode_final_parameters[n, 21],
      a42 = NativeCommunities_stode_final_parameters[n, 22],
      a43 = NativeCommunities_stode_final_parameters[n, 23],
      a44 = NativeCommunities_stode_final_parameters[n, 24],
      a45 = NativeCommunities_stode_final_parameters[n, 25],
      a51 = NativeCommunities_stode_final_parameters[n, 26],
      a52 = NativeCommunities_stode_final_parameters[n, 27],
      a53 = NativeCommunities_stode_final_parameters[n, 28],
      a54 = NativeCommunities_stode_final_parameters[n, 29],
      a55 = NativeCommunities_stode_final_parameters[n, 30],
      muX6 = SKparamsForSKreassess[n, 31],
      a16 = SKparamsForSKreassess[n, 32],
      a61 = SKparamsForSKreassess[n, 33],
      a26 = SKparamsForSKreassess[n, 34],
      a62 = SKparamsForSKreassess[n, 35],
      a36 = SKparamsForSKreassess[n, 36],
      a63 = SKparamsForSKreassess[n, 37],
      a46 = SKparamsForSKreassess[n, 38],
      a64 = SKparamsForSKreassess[n, 39],
      a56 = SKparamsForSKreassess[n, 40],
      a65 = SKparamsForSKreassess[n, 41],
      a66 = SKparamsForSKreassess[n, 42])
    
    state <- c(X1 = NativeCommunities_stode_final_states[n, 1],
               X2 = NativeCommunities_stode_final_states[n, 2],
               X3 = NativeCommunities_stode_final_states[n, 3],
               X4 = NativeCommunities_stode_final_states[n, 4],
               X5 = NativeCommunities_stode_final_states[n, 5],
               X6 = start_amt)
    
    # the non-linear solver
    out_4 <- stode(y = state, fun = NonLinearGLV6, parms = parametersSingleKiller, pos = TRUE)
    
    # calculate jacobian for eigenvalue evaluation of stability
    jac <- jacobian.full(y = c(out_4$y[1], out_4$y[2], out_4$y[3], out_4$y[4], out_4$y[5], out_4$y[6]), func = NonLinearGLV6, parms = parametersSingleKiller)
    eig <- eigen(jac)
    
    # fills matrices from each nat com simulation
    final_states_singlekiller[n, ] = out_4$y
    final_parameters_singlekiller[n, ] <- parametersSingleKiller
    final_eigs_singlekiller[n, ] <- Re(eig$values)
    
    n <- n + 1
  }

  # create directories for storing newly evolved killer comm data
  dir.create(paste0(wrkdir, "evolution_of_PSC_", PSC_number, "/Round_", z, "_data/KillerComm_stode_community_eigs"))
  dir.create(paste0(wrkdir, "evolution_of_PSC_", PSC_number, "/Round_", z, "_data/KillerComm_stode_community_states"))
  dir.create(paste0(wrkdir, "evolution_of_PSC_", PSC_number, "/Round_", z, "_data/KillerComm_stode_community_parameters"))
  
  #writing single killer parameters file into evolution round directory
  csvFileName <- paste("SingleKiller_stode_final_parameters.csv")
  csvFilePathSKparams <- paste0(wrkdir, "evolution_of_PSC_", PSC_number, "/Round_", z, "_data/", csvFileName) 
  write.table(final_parameters_singlekiller, file = csvFilePathSKparams, sep = ",", col.names = T, row.names = T)
  
  #writing single killer states file into evolution round directory
  csvFileName <- paste("SingleKiller_stode_final_states.csv")
  csvFilePathSKstates <-paste0(wrkdir, "evolution_of_PSC_", PSC_number, "/Round_", z, "_data/", csvFileName) 
  write.table(final_states_singlekiller, file = csvFilePathSKstates, sep = ",", col.names = F, row.names = F)
  
  #writing single killer eigs file into evolution round directory
  csvFileName <- paste("SingleKiller_stode_final_eigs.csv")
  csvFilePathSKeigs <- paste0(wrkdir, "evolution_of_PSC_", PSC_number, "/Round_", z, "_data/", csvFileName)
  write.table(final_eigs_singlekiller, file = csvFilePathSKeigs, sep = ",", col.names = F, row.names = F)
  
  # rename. should fix this for consistency
  SingleKiller_stode_final_parameters <- final_parameters_singlekiller
  
  #begin KillerComm community generation, reset loop counters
  n <- 1  # loop counter for number of native + killer combos
  p <- 1  # loop counter for number of KillerComms
  
  print("evolved PSC generation...")
  
  # set seed to ensure a new set of "random" numbers
  set.seed(topPSC + z)
  
  # begins loop that will be repeated for each killercomm
  for (p in 1:num_of_KillerComm){
    
    # create killer comm parameter matrix for later export
    final_parameters_KillerComm <- matrix(ncol = (num_of_species*num_of_species+num_of_species), nrow = num_of_comms)
    colnames(final_parameters_KillerComm) <- c("muX1",
                                               "muX2",
                                               "muX3",
                                               "muX4",
                                               "muX5",
                                               "a11",
                                               "a12",
                                               "a13",
                                               "a14",
                                               "a15",
                                               "a21",
                                               "a22",
                                               "a23",
                                               "a24",
                                               "a25",
                                               "a31",
                                               "a32",
                                               "a33",
                                               "a34",
                                               "a35",
                                               "a41",
                                               "a42",
                                               "a43",
                                               "a44",
                                               "a45",
                                               "a51",
                                               "a52",
                                               "a53",
                                               "a54",
                                               "a55",
                                               "muX6",
                                               "a16",
                                               "a61",
                                               "a26",
                                               "a62",
                                               "a36",
                                               "a63",
                                               "a46",
                                               "a64",
                                               "a56",
                                               "a65",
                                               "a66",
                                               "muX7",
                                               "muX8",
                                               "muX9",
                                               "muX10",
                                               "muX11",
                                               "a67",
                                               "a68",
                                               "a69",
                                               "a610",
                                               "a611",
                                               "a116",
                                               "a106",
                                               "a96",
                                               "a86",
                                               "a76",
                                               "a17",
                                               "a18",
                                               "a19",
                                               "a110",
                                               "a1_11",
                                               "a27",
                                               "a28",
                                               "a29",
                                               "a210",
                                               "a211",
                                               "a37",
                                               "a38",
                                               "a39",
                                               "a310",
                                               "a311",
                                               "a47",
                                               "a48",
                                               "a49",
                                               "a410",
                                               "a411",
                                               "a57",
                                               "a58",
                                               "a59",
                                               "a510",
                                               "a511",
                                               "a71",
                                               "a72",
                                               "a73",
                                               "a74",
                                               "a75",
                                               "a77",
                                               "a78",
                                               "a79",
                                               "a710",
                                               "a711",
                                               "a81",
                                               "a82",
                                               "a83",
                                               "a84",
                                               "a85",
                                               "a87",
                                               "a88",
                                               "a89",
                                               "a810",
                                               "a811",
                                               "a91",
                                               "a92",
                                               "a93",
                                               "a94",
                                               "a95",
                                               "a97",
                                               "a98",
                                               "a99",
                                               "a910",
                                               "a911",
                                               "a101",
                                               "a102",
                                               "a103",
                                               "a104",
                                               "a105",
                                               "a107",
                                               "a108",
                                               "a109",
                                               "a1010",
                                               "a1011",
                                               "a11_1",
                                               "a112",
                                               "a113",
                                               "a114",
                                               "a115",
                                               "a117",
                                               "a118",
                                               "a119",
                                               "a1110",
                                               "a1111")
    # create killer comm final states matrix for later export
    final_states_KillerComm <- matrix(ncol = num_of_species, nrow = num_of_comms)
    
    # create killer comm eigs matrix for later export
    final_eigs_KillerComm <- matrix(ncol=num_of_species, nrow=num_of_comms)
    
    # generate new set of parameters that will be replaced by important parameters identified above
    KCparms <- c(
      #KillerComm growth rates
      muX7 = abs(rnorm(1, mean = mu_mean, sd = mu_sd)),
      muX8 = abs(rnorm(1, mean = mu_mean, sd = mu_sd)),
      muX9 = abs(rnorm(1, mean = mu_mean, sd = mu_sd)),
      muX10 = abs(rnorm(1, mean = mu_mean, sd = mu_sd)),
      muX11 = abs(rnorm(1, mean = mu_mean, sd = mu_sd)),
      #killerComm interaction parameters with Killer
      a67 = rnorm(1, mean = param_mean, sd = param_sd),
      a68 = rnorm(1, mean = param_mean, sd = param_sd),
      a69 = rnorm(1, mean = param_mean, sd = param_sd),
      a610 = rnorm(1, mean = param_mean, sd = param_sd),
      a611 = rnorm(1, mean = param_mean, sd = param_sd),
      a116 = rnorm(1, mean = param_mean, sd = param_sd),
      a106 = rnorm(1, mean = param_mean, sd = param_sd),
      a96 = rnorm(1, mean = param_mean, sd = param_sd),
      a86 = rnorm(1, mean = param_mean, sd = param_sd),
      a76 = rnorm(1, mean = param_mean, sd = param_sd),
      # NativeComm Interactions with Killer Comm
      a17 = rnorm(1, mean = param_mean, sd = param_sd),
      a18 = rnorm(1, mean = param_mean, sd = param_sd),
      a19 = rnorm(1, mean = param_mean, sd = param_sd),
      a110 = rnorm(1, mean = param_mean, sd = param_sd),
      a1_11 = rnorm(1, mean = param_mean, sd = param_sd),
      a27 = rnorm(1, mean = param_mean, sd = param_sd),
      a28 = rnorm(1, mean = param_mean, sd = param_sd),
      a29 = rnorm(1, mean = param_mean, sd = param_sd),
      a210 = rnorm(1, mean = param_mean, sd = param_sd),
      a211 = rnorm(1, mean = param_mean, sd = param_sd),
      a37 = rnorm(1, mean = param_mean, sd = param_sd),
      a38 = rnorm(1, mean = param_mean, sd = param_sd),
      a39 = rnorm(1, mean = param_mean, sd = param_sd),
      a310 = rnorm(1, mean = param_mean, sd = param_sd),
      a311 = rnorm(1, mean = param_mean, sd = param_sd),
      a47 = rnorm(1, mean = param_mean, sd = param_sd),
      a48 = rnorm(1, mean = param_mean, sd = param_sd),
      a49 = rnorm(1, mean = param_mean, sd = param_sd),
      a410 = rnorm(1, mean = param_mean, sd = param_sd),
      a411 = rnorm(1, mean = param_mean, sd = param_sd),
      a57 = rnorm(1, mean = param_mean, sd = param_sd),
      a58 = rnorm(1, mean = param_mean, sd = param_sd),
      a59 = rnorm(1, mean = param_mean, sd = param_sd),
      a510 = rnorm(1, mean = param_mean, sd = param_sd),
      a511 = rnorm(1, mean = param_mean, sd = param_sd),
      #KillerComm interactions with native community
      a71 = rnorm(1, mean = param_mean, sd = param_sd),
      a72 = rnorm(1, mean = param_mean, sd = param_sd),
      a73 = rnorm(1, mean = param_mean, sd = param_sd),
      a74 = rnorm(1, mean = param_mean, sd = param_sd),
      a75 = rnorm(1, mean = param_mean, sd = param_sd),
      a77 = abs(rnorm(1, mean = param_mean, sd = param_sd))*-1,
      a78 = rnorm(1, mean = param_mean, sd = param_sd),
      a79 = rnorm(1, mean = param_mean, sd = param_sd),
      a710 = rnorm(1, mean = param_mean, sd = param_sd),
      a711 = rnorm(1, mean = param_mean, sd = param_sd),
      a81 = rnorm(1, mean = param_mean, sd = param_sd),
      a82 = rnorm(1, mean = param_mean, sd = param_sd),
      a83 = rnorm(1, mean = param_mean, sd = param_sd),
      a84 = rnorm(1, mean = param_mean, sd = param_sd),
      a85 = rnorm(1, mean = param_mean, sd = param_sd),
      a87 = rnorm(1, mean = param_mean, sd = param_sd),
      a88 = abs(rnorm(1, mean = param_mean, sd = param_sd))*-1,
      a89 = rnorm(1, mean = param_mean, sd = param_sd),
      a810 = rnorm(1, mean = param_mean, sd = param_sd),
      a811 = rnorm(1, mean = param_mean, sd = param_sd),
      a91 = rnorm(1, mean = param_mean, sd = param_sd),
      a92 = rnorm(1, mean = param_mean, sd = param_sd),
      a93 = rnorm(1, mean = param_mean, sd = param_sd),
      a94 = rnorm(1, mean = param_mean, sd = param_sd),
      a95 = rnorm(1, mean = param_mean, sd = param_sd),
      a97 = rnorm(1, mean = param_mean, sd = param_sd),
      a98 = rnorm(1, mean = param_mean, sd = param_sd),
      a99 = abs(rnorm(1, mean = param_mean, sd = param_sd))*-1,
      a910 = rnorm(1, mean = param_mean, sd = param_sd),
      a911 = rnorm(1, mean = param_mean, sd = param_sd),
      a101 = rnorm(1, mean = param_mean, sd = param_sd),
      a102 = rnorm(1, mean = param_mean, sd = param_sd),
      a103 = rnorm(1, mean = param_mean, sd = param_sd),
      a104 = rnorm(1, mean = param_mean, sd = param_sd),
      a105 = rnorm(1, mean = param_mean, sd = param_sd),
      a107 = rnorm(1, mean = param_mean, sd = param_sd),
      a108 = rnorm(1, mean = param_mean, sd = param_sd),
      a109 = rnorm(1, mean = param_mean, sd = param_sd),
      a1010 = abs(rnorm(1, mean = param_mean, sd = param_sd))*-1,
      a1011 = rnorm(1, mean = param_mean, sd = param_sd),
      a11_1 = rnorm(1, mean = param_mean, sd = param_sd),
      a112 = rnorm(1, mean = param_mean, sd = param_sd),
      a113 = rnorm(1, mean = param_mean, sd = param_sd),
      a114 = rnorm(1, mean = param_mean, sd = param_sd),
      a115 = rnorm(1, mean = param_mean, sd = param_sd),
      a117 = rnorm(1, mean = param_mean, sd = param_sd),
      a118 = rnorm(1, mean = param_mean, sd = param_sd),
      a119 = rnorm(1, mean = param_mean, sd = param_sd),
      a1110 = rnorm(1, mean = param_mean, sd = param_sd),
      a1111 = abs(rnorm(1, mean = param_mean, sd = param_sd))*-1)
    
    # loop to replace fixed values NOT WORKING
    s <- 1
    for (s in 1:nrow(importantParms)) {
      KCparms <- replace(KCparms, as.character(importantParms[s, 7]), as.numeric(importantParms[s, 2]))
      s <- s + 1
    }
    
    # begins loop that will be repeated for each combo of killer + native community
    for (n in 1:num_of_comms) {
      
      # reading in parameters from single killer model and native communities
      NCandSKparms <- c(
        muX1 = NativeCommunities_stode_final_parameters[n, 1], 
        muX2 = NativeCommunities_stode_final_parameters[n, 2], 
        muX3 = NativeCommunities_stode_final_parameters[n, 3], 
        muX4 = NativeCommunities_stode_final_parameters[n, 4], 
        muX5 = NativeCommunities_stode_final_parameters[n, 5], 
        a11 = NativeCommunities_stode_final_parameters[n, 6], 
        a12 = NativeCommunities_stode_final_parameters[n, 7], 
        a13 = NativeCommunities_stode_final_parameters[n, 8], 
        a14 = NativeCommunities_stode_final_parameters[n, 9], 
        a15 = NativeCommunities_stode_final_parameters[n, 10], 
        a21 = NativeCommunities_stode_final_parameters[n, 11], 
        a22 = NativeCommunities_stode_final_parameters[n, 12], 
        a23 = NativeCommunities_stode_final_parameters[n, 13], 
        a24 = NativeCommunities_stode_final_parameters[n, 14], 
        a25 = NativeCommunities_stode_final_parameters[n, 15], 
        a31 = NativeCommunities_stode_final_parameters[n, 16], 
        a32 = NativeCommunities_stode_final_parameters[n, 17], 
        a33 = NativeCommunities_stode_final_parameters[n, 18], 
        a34 = NativeCommunities_stode_final_parameters[n, 19], 
        a35 = NativeCommunities_stode_final_parameters[n, 20], 
        a41 = NativeCommunities_stode_final_parameters[n, 21], 
        a42 = NativeCommunities_stode_final_parameters[n, 22], 
        a43 = NativeCommunities_stode_final_parameters[n, 23], 
        a44 = NativeCommunities_stode_final_parameters[n, 24], 
        a45 = NativeCommunities_stode_final_parameters[n, 25], 
        a51 = NativeCommunities_stode_final_parameters[n, 26], 
        a52 = NativeCommunities_stode_final_parameters[n, 27], 
        a53 = NativeCommunities_stode_final_parameters[n, 28], 
        a54 = NativeCommunities_stode_final_parameters[n, 29], 
        a55 = NativeCommunities_stode_final_parameters[n, 30], 
        #single killer parameters 
        muX6 = as.numeric(SingleKiller_stode_final_parameters[n, 31]), 
        a16 = as.numeric(SingleKiller_stode_final_parameters[n, 32]), 
        a61 = as.numeric(SingleKiller_stode_final_parameters[n, 33]), 
        a26 = as.numeric(SingleKiller_stode_final_parameters[n, 34]), 
        a62 = as.numeric(SingleKiller_stode_final_parameters[n, 35]), 
        a36 = as.numeric(SingleKiller_stode_final_parameters[n, 36]), 
        a63 = as.numeric(SingleKiller_stode_final_parameters[n, 37]), 
        a46 = as.numeric(SingleKiller_stode_final_parameters[n, 38]), 
        a64 = as.numeric(SingleKiller_stode_final_parameters[n, 39]), 
        a56 = as.numeric(SingleKiller_stode_final_parameters[n, 40]), 
        a65 = as.numeric(SingleKiller_stode_final_parameters[n, 41]), 
        a66 = as.numeric(SingleKiller_stode_final_parameters[n, 42]))
      
      # sets starting abundance to the final states obtained from the previous models and adds killercomm abundances
      state <- c(X1 = NativeCommunities_stode_final_states[n, 1],
                 X2 = NativeCommunities_stode_final_states[n, 2],
                 X3 = NativeCommunities_stode_final_states[n, 3],
                 X4 = NativeCommunities_stode_final_states[n, 4],
                 X5 = NativeCommunities_stode_final_states[n, 5],
                 X6 = start_amt,
                 X7 = start_amt,
                 X8 = start_amt,
                 X9 = start_amt,
                 X10 = start_amt,
                 X11 = start_amt)

      # assigns growth rates for all species and interaction coefficients as parameters for nonlinear solver
      parametersKillerComm <- c(NCandSKparms, KCparms)
     
      # the non-linear solver
      out_4 <- stode(y = state, fun = NonLinearGLV11, parms = parametersKillerComm, pos = TRUE)
      
      # eigen-thingy for stability
      jac <- jacobian.full(y = c(out_4$y[1], out_4$y[2], out_4$y[3], out_4$y[4], out_4$y[5], out_4$y[6], out_4$y[7], out_4$y[8], out_4$y[9], out_4$y[10], out_4$y[11]), func = NonLinearGLV11, parms = parametersKillerComm)
      eig <- eigen(jac)
      
      # fills the matrices
      final_states_KillerComm[n, ] = out_4$y
      final_parameters_KillerComm[n, ] <- parametersKillerComm
      final_eigs_KillerComm[n, ] <- Re(eig$values)
      
      n = n + 1
    }
    
    # write the killer comm parameter files
    csvFileName <- paste("KillerComm_",p,"_stode_final_parameters.csv",sep="")
    csvFilePath <- paste0(wrkdir, "evolution_of_PSC_", PSC_number, "/Round_", z, "_data/KillerComm_stode_community_parameters/",csvFileName) 
    write.table(final_parameters_KillerComm, file = csvFilePath, sep = ",", col.names = T, row.names = F)
    
    # write the killer comm final states files
    csvFileName2 <- paste("KillerComm_",p,"_stode_final_states.csv",sep="")
    csvFilePath2 <- paste0(wrkdir, "evolution_of_PSC_", PSC_number, "/Round_", z, "_data/KillerComm_stode_community_states/",csvFileName2) 
    write.table(final_states_KillerComm, file = csvFilePath2, sep = ",", col.names = F, row.names = F)
    
    # write the killer comm eigs files
    csvFileName3 <- paste("KillerComm_",p,"_stode_final_eigs.csv",sep="")
    csvFilePath3 <- paste0(wrkdir, "evolution_of_PSC_", PSC_number, "/Round_", z, "_data/KillerComm_stode_community_eigs/",csvFileName3) 
    write.table(final_eigs_KillerComm, file = csvFilePath3, sep = ",", col.names = F, row.names = F)
    
    if (p %% 10 == 0){
      print(p)
    }
    
    p = p + 1
    n = 1
  }
  
  #end PSC generation, begin PSC analysis
  
  KillerCommMetrics <- matrix(nrow = num_of_KillerComm, ncol = 22)
  col_names_KillerCommMetrics <- rbind("Average of Total Absolute Community Abundance across NCs","Average of Probiotic Support Community Absolute Community Abundance across NCs","Average of Target Species Absolute Community Abundance across NCs","Average of NC Absolute Community Abundance across NCs","Average Probiotic Support Community Relative Community Abundance","Average of Target Species Relative Community Abundance","Average of NC Combined Relative Community Abundance","# of NC Target Species Diff > 0 (total abundance basis)","# of NC Target Species Diff < 0 (total abundance basis)","# of NC Target Species Diff > 0 (relative abundance basis)","# of NC Target Species Diff < 0 (relative abundance basis)","# of NC Target Species Total Abundance > 0.001","# of NC Target Species Total Abundance > 0.01","# of NC Target Species Total Abundance > 0.1","# of NC Target Species Relative Abundance > 0.001","# of NC Target Species Relative Abundance > 0.01","# of NC Target Species Relative Abundance > 0.1","Average mean evenness of final communities across all NCs","# of NC with all (-) eigenvalues","Probiotic Support Community Number","# of NC whereby target species abundance = 0 in single target speices inoculation AND target species abundance > 0 in PSC inoculation","# of NC whereby target species abundance > 0 in single target speices inoculation AND target species abundance > in PSC inoculation than in single target species inoculation")
  colnames(KillerCommMetrics) <- col_names_KillerCommMetrics
  
  SingleKillerMetrics <- matrix(ncol = 13, nrow = 1)
  
  #importing Single Killer Final States and Eigs
  SingleKillerEigs<- final_eigs_singlekiller #read.csv(paste("~/Documents/SynCom_Modeling/",file_number,"/SingleKiller_stode_final_eigs.csv",sep=""), header=FALSE)
  
  SingleKillerFinStates <- final_states_singlekiller #read.csv(paste("~/Documents/SynCom_Modeling/",file_number,"/SingleKiller_stode_final_states.csv",sep=""), header=FALSE)
  sumTotalAbund <- rowSums(SingleKillerFinStates, na.rm = TRUE) #total single killer community abundance
  sumNativeAbund <- rowSums(SingleKillerFinStates[, 1:5], na.rm = TRUE) #native community abundance from 
  SingleKillerFinStates <- cbind(SingleKillerFinStates, sumTotalAbund, sumNativeAbund)
  
  SingleKillerRelStates <- matrix(ncol = 18, nrow = num_of_comms)
  SingleKillerRelStates[, 1] = SingleKillerFinStates[, 1] / SingleKillerFinStates[, 7]
  SingleKillerRelStates[, 2] = SingleKillerFinStates[, 2] / SingleKillerFinStates[, 7]
  SingleKillerRelStates[, 3] = SingleKillerFinStates[, 3] / SingleKillerFinStates[, 7]
  SingleKillerRelStates[, 4] = SingleKillerFinStates[, 4] / SingleKillerFinStates[, 7]
  SingleKillerRelStates[, 5] = SingleKillerFinStates[, 5] / SingleKillerFinStates[, 7]
  SingleKillerRelStates[, 6] = SingleKillerFinStates[, 6] / SingleKillerFinStates[, 7]
  SingleKillerMetrics[, 1] <- sum(SingleKillerFinStates[, 6] > 0.001, na.rm = TRUE) # of native comms target species total abundance > 0.001
  SingleKillerMetrics[, 2] <- sum(SingleKillerFinStates[, 6] > 0.01, na.rm = TRUE) # of native comms target species total abundance > 0.01
  SingleKillerMetrics[, 3] <- sum(SingleKillerFinStates[, 6] > 0.1, na.rm = TRUE) # of native comms target species total abundance > 0.1
  SingleKillerMetrics[, 4] <- sum(SingleKillerRelStates[, 6] > 0.001, na.rm = TRUE) # of native comms target species total abundance > 0.001
  SingleKillerMetrics[, 5] <- sum(SingleKillerRelStates[, 6] > 0.01, na.rm = TRUE) # of native comms target species total abundance > 0.01
  SingleKillerMetrics[, 6] <- sum(SingleKillerRelStates[, 6] > 0.1, na.rm = TRUE) # of native comms target species total abundance > 0.1
  
  #calculating evenness for each Single Killer native comm (p*ln(p))
  n <- 1
  for(n in 1:6) {
    SingleKillerRelStates[, 6 + n] <- SingleKillerRelStates[, n] * (log(SingleKillerRelStates[, n]))
  }
  
  SingleKillerRelStates[, 13] <- rowSums(SingleKillerRelStates[, 7:12], na.rm = TRUE) #summing p*ln(p) values across all 11 species
  SingleKillerRelStates[, 14] <- (SingleKillerRelStates[, 13]) * -1 #multiply sum by -1
  SingleKillerRelStates[, 15] <- SingleKillerRelStates[, 14] / log(6) # divide by log(6) - column 15 is evenness for each community
  SingleKillerRelStates[, 16] <- SingleKillerRelStates[, 6] #rel abundance target species single killer
  SingleKillerRelStates[, 17] <- rowSums(SingleKillerRelStates[, 1:5], na.rm = TRUE) #rel abundance native community single killer
  
  SingleKillerMetrics[, 7] <- mean(SingleKillerRelStates[, 15], na.rm = TRUE) #mean evenness of single killer final equilibrated communities
  SingleKillerMetrics[, 8] <- length(which(rowSums(SingleKillerEigs < 0) == 6))  #number of single killer communities with all negative eigenvalues
  SingleKillerMetrics[, 9] <- mean(SingleKillerFinStates[, 6], na.rm = TRUE) #single killer target species total abundance
  SingleKillerMetrics[, 10] <- mean(SingleKillerFinStates[, 7], na.rm = TRUE) #single killer total community total abundance
  SingleKillerMetrics[, 11] <- mean(SingleKillerFinStates[, 8], na.rm = TRUE) #single killer native community total abundance
  SingleKillerMetrics[, 12] <- mean(SingleKillerRelStates[, 6], na.rm = TRUE) #single killer target species relative abuundance
  SingleKillerMetrics[, 13] <- mean(SingleKillerRelStates[, 17], na.rm = TRUE) #single killer native community relative abundance
  
  p <- 1
  for(p in 1:num_of_KillerComm) {
    #read in final states for KillerComm = p
    csvFileName <- paste0("KillerComm_",p,"_stode_final_states.csv")
    csvFilePath <- paste0(wrkdir, "evolution_of_PSC_", PSC_number, "/Round_", z, "_data/KillerComm_stode_community_states/",csvFileName) 
    commStates<-read.table(file = csvFilePath, sep=",")
    #read in final eigenvalues for KillerComm = p
    csvFileName2 <- paste0("KillerComm_",p,"_stode_final_eigs.csv")
    csvFilePath2 <- paste0(wrkdir, "evolution_of_PSC_", PSC_number, "/Round_", z, "_data/KillerComm_stode_community_eigs/",csvFileName2) 
    commEigs<-read.table(file = csvFilePath2, sep=",")
    
    #make zero vectors to tack onto as columns to commStates
    if(num_of_species!=11){
      num_of_columns_to_add<-11-num_of_species
      zeroColumns<-matrix(0L, ncol=num_of_columns_to_add, nrow=num_of_comms)
      commStates<-cbind(commStates, zeroColumns)
    }
    
    #Begin community metric calculations
    sumTotalAbund <- rowSums(commStates, na.rm = TRUE) #summing total absolute abundance of each final community
    sumKillerCommTotalAbund <- rowSums(commStates[, 7:11], na.rm = TRUE)
    targetTotalAbund <- commStates[, 6]
    nativeCommTotalAbund <- rowSums(commStates[, 1:5], na.rm = TRUE)
    KillerCommTotalDiffs <- commStates[, 6] - SingleKillerFinStates[, 6]
    commStates <- cbind(commStates, sumTotalAbund, sumKillerCommTotalAbund, targetTotalAbund, nativeCommTotalAbund, KillerCommTotalDiffs)
    KillerCommMetrics[p, 1] <- mean(commStates[, 12]) #average of Total Absolute Community Abundance into KillerCommMetrics
    KillerCommMetrics[p, 2] <- mean(commStates[, 13]) #average of KillerComm Absolute Community Abundance into KillerCommMetrics
    KillerCommMetrics[p, 3] <- mean(commStates[, 14]) #average of Target Species Absolute Community Abundance into KillerCommMetrics
    KillerCommMetrics[p, 4] <- mean(commStates[, 15]) #average of NativeComm Absolute Community Abundance into KillerCommMetrics
    
    #calculating relative abundances of each community member
    commStatesRel <- matrix(ncol = 25, nrow = num_of_comms)
    commStatesRel[, 1] <- commStates[, 1] / commStates[, 12]
    commStatesRel[, 2] <- commStates[, 2] / commStates[, 12]
    commStatesRel[, 3] <- commStates[, 3] / commStates[, 12]
    commStatesRel[, 4] <- commStates[, 4] / commStates[, 12]
    commStatesRel[, 5] <- commStates[, 5] / commStates[, 12]
    commStatesRel[, 6] <- commStates[, 6] / commStates[, 12]
    commStatesRel[, 7] <- commStates[, 7] / commStates[, 12]
    commStatesRel[, 8] <- commStates[, 8] / commStates[, 12]
    commStatesRel[, 9] <- commStates[, 9] / commStates[, 12]
    commStatesRel[, 10] <- commStates[, 10] / commStates[, 12]
    commStatesRel[, 11] <- commStates[, 11] / commStates[, 12]
    sumKillerCommRelAbund <- rowSums(commStatesRel[, 7:11])
    targetRelAbund <- commStatesRel[, 6]
    nativeCommRelAbund <- rowSums(commStatesRel[, 1:5])
    KillerCommRelDiffs <- commStatesRel[, 6] - SingleKillerRelStates[, 6]
    #calculating evenness for each final comm
    
    n <- 1
    #calculating p*ln(p) for each species
    for(n in 1:num_of_species) {
      commStatesRel[, 11 + n] <- commStatesRel[, n] * (log(commStatesRel[, n]))
    } 
    
    commStatesRel[, 23] <- rowSums(commStatesRel[, 12:22], na.rm = TRUE) #summing p*ln(p) values across all species 
    commStatesRel[, 24] <- (commStatesRel[, 23]) * -1 #multiply sum by -1
    commStatesRel[, 25] <- commStatesRel[, 24] / log(num_of_species) # divide by log(11) - this is final comm evenness 
    
    #entering metrics into KillerCommMetrics
    commStatesRel <- cbind(commStatesRel, sumKillerCommRelAbund, targetRelAbund, nativeCommRelAbund, KillerCommRelDiffs)
    commStatesRel[is.na(commStatesRel)] <- 0
    KillerCommMetrics[p, 5] <- mean(commStatesRel[, 26], na.rm = TRUE) #average of KillerComm Combined Relative Community Abundance into KillerCommMetrics
    KillerCommMetrics[p, 6] <- mean(commStatesRel[, 27], na.rm = TRUE) #average of target species Relative Community Abundance into KillerCommMetrics
    KillerCommMetrics[p, 7] <- mean(commStatesRel[, 28], na.rm = TRUE) #average of NativeComm Combined Relative Community Abundance into KillerCommMetrics
    KillerCommMetrics[p, 8] <- sum(commStates[, 16] > 0, na.rm = TRUE)# of Native Comms KillerComm target species Diff > 0 (total)
    KillerCommMetrics[p, 9] <- sum(commStates[, 16] < 0, na.rm = TRUE)# of Native Comms KillerComm target species Diff < 0 (total)
    KillerCommMetrics[p, 10] <- sum(commStatesRel[, 29] > 0, na.rm = TRUE) # of Native Comms KillerComm target species Diff > 0 (relative)
    KillerCommMetrics[p, 11] <- sum(commStatesRel[, 29] < 0, na.rm = TRUE) # of Native Comms KillerComm target species Diff < 0 (relative)
    KillerCommMetrics[p, 12] <- sum(commStates[, 6] > 0.001, na.rm = TRUE) # of native comms target species total abundance > 0.001
    KillerCommMetrics[p, 13] <- sum(commStates[, 6] > 0.01, na.rm = TRUE) # of native comms target species total abundance > 0.01
    KillerCommMetrics[p, 14] <- sum(commStates[, 6] > 0.1, na.rm = TRUE) # of native comms target species total abundance > 0.1
    KillerCommMetrics[p, 15] <- sum(commStatesRel[, 6] > 0.001, na.rm = TRUE) # of comms target species relative abundance > 0.001
    KillerCommMetrics[p, 16] <- sum(commStatesRel[, 6] > 0.01, na.rm = TRUE) # of comms target species relative abundance > 0.01
    KillerCommMetrics[p, 17] <- sum(commStatesRel[, 6] > 0.1, na.rm = TRUE) # of comms target species relative abundance > 0.1
    KillerCommMetrics[p, 18] <- mean(commStatesRel[, 25], na.rm = TRUE) #mean evenness of final community
    num_of_neg_eigs <- rowSums(commEigs < 0)
    KillerCommMetrics[p, 19] <- length(which(num_of_neg_eigs == num_of_species)) # number of communities with all negative eigenvalues
    KillerCommMetrics[p, 20] <- p #records probiotic support community number for row 
    
    # num of NC [0 in single species inoc] -> [>0 in probiotic support community inoc]
    #1. [Single Killer Fin States - species 6] AND [final States Killer Comm - species 6] into same matrix, two columns
    perform_metrics_table <- matrix(nrow = num_of_comms, ncol = 5)
    perform_metrics_table[, 1] <- commStates[, 6]
    perform_metrics_table[, 2] <- SingleKillerFinStates[, 6]
    #2. column 3, give a 1 if single killer Fin state = 0
    perform_metrics_table[, 3] <- ifelse(perform_metrics_table[, 2] == 0, 1, 0)
    #3. column 4, give a 1 if [killercomm-singleinoc]>0
    perform_metrics_table[, 4] <- ifelse(perform_metrics_table[, 1] > perform_metrics_table[, 2], 1, 0)
    #4. column 5, add two columns
    perform_metrics_table[, 5] <- perform_metrics_table[, 3] + perform_metrics_table[, 4]
    #5. NC of [0] -> [>0] is sums(column 5 = 2) - number of rows whose column #5 is = 2
    NCs_zero_to_pos <- sum(perform_metrics_table[, 5] == 2)
    #6. Put this # of NC inot KillerCommMetrics[p,21]
    KillerCommMetrics[p, 21] <- NCs_zero_to_pos # of NC whereby target species abundance = 0 in single target speices inoculation AND target species abundance > 0 in PSC inoculation
    
    perform_metrics_table2 <- matrix(nrow = num_of_comms, ncol = 5)
    perform_metrics_table2[, 1] <- commStates[, 6] #final states from PSC
    perform_metrics_table2[, 2] <- SingleKillerFinStates[, 6] #final states from single target species inoculation
    # of NC [>0 in single species inoc] -> [>>0 in probiotic support community inoc] 
    #1. [Single Killer Fin States - species 6] AND [final States Killer Comm - species 6] into same matrix, two columns
    #2. column 3, give a 1 if single killer Fin state > 0
    perform_metrics_table2[, 3] <- ifelse(perform_metrics_table2[, 2] > 0, 1, 0)
    #3. column 4, give a 1 if [killercomm-singleinoc]>0
    perform_metrics_table2[, 4] <- ifelse(perform_metrics_table2[, 1] > perform_metrics_table2[, 2], 1, 0)
    #4. column 5, add two columns
    perform_metrics_table2[, 5] <- perform_metrics_table2[, 3] + perform_metrics_table2[, 4]
    #5. NC of [0] -> [>0] is sums(column 5 = 2) - number of rows whose column #5 is = 2
    NCs_pos_to_greater <- sum(perform_metrics_table2[, 5] == 2)
    #6. Put this # of NC inot ParameterSensitivity[p,22]
    KillerCommMetrics[p, 22] <- NCs_pos_to_greater # of NC whereby target species abundance > 0 in single target speices inoculation AND target species abundance > in PSC inoculation than in single target species inoculation
    
    rm(perform_metrics_table)
    rm(perform_metrics_table2)
    rm(commStatesRel)
    
    p = p + 1
  }
  
  averageMetrics <- colMeans(KillerCommMetrics, na.rm = TRUE)
  maxMetrics <- apply(KillerCommMetrics, 2, max, na.rm = TRUE) #Finds max of each column in KillerCommMetrics matrix
  minMetrics <- apply(KillerCommMetrics, 2, min, na.rm = TRUE) #Finds min of each column in KillerCommMetrics matrix
  
  #combining metrics tables, generating 'combinedMetrics' and writing to file directory for set
  combinedMetrics <- matrix(nrow = 1, ncol = 77)
  col_names_combinedMetrics1 <- rbind("file_set number","Average Total Abundance – Whole Community","Average Total Abundance – Probiotic Support Community","Average Total Abundance – Target Species","Average Total Abundance – Native Community","AverageRelative Abundance – Probiotic Support Community","Average Relative Abundance – Target Species","Average Relative Abundance – Native Community","Average # of NC Target Species Diff > 0 (total abundance, average of all probiotic support communities)","Average # of NC Target Species Diff < 0 (total abundance, average of all probiotic support communities)","Average # of NC Target Species Diff > 0 (relative abundance, average of all probiotic support communities)","Average # of NC Target Species Diff < 0 (relative abundance, average of all probiotic support communities)","Average # of NC target species total abundance > 0.001","Average # of NC target species total abundance > 0.01","Average # of NC target species total abundance > 0.1","Average # of NC target species relative abundance > 0.001","Average # of NC target species relative abundance > 0.01","Average # of NC target species relative abundance > 0.1","Average of Average Community Evenness","Average # of NC with all (-) eigenvalues","Average # of NCs 0 -> >0","Average # of NCs >0 -> >>0")
  col_names_combinedMetrics2 <- rbind("Maximum Total Abundance - Whole Community","Maximum Total Abundance - Probiotic Support Community","Maximum Total Abundance - Target Species","Maximum Total Abundance - Native Community","Maximum Relative Abundance - Probiotic Support Community","Maximum Relative Abundance - Target Species","Maximum Relative Abundance - Native Community","Maximum # of NC Target Species Diff > 0 (total abundance basis)","Maximum # of NC Target Species Diff < 0 (total abundance basis)","Maximum # of NC Target Species Diff > 0 (relative abundance basis)","Maximum # of NC Target Species Diff < 0 (relative abundance basis)","Maximum # of NC target species total abundance > 0.001","Maximum # of NC target species total abundance > 0.01","Maximum # of NC target species total abundance > 0.1","Maximum # of NC target species relative abundance > 0.001","Maximum # of NC target species relative abundance > 0.01","Maximum # of NC target species relative abundance > 0.1","Maxiumum Average Final Combined Community Evenness","Maximum # of NC with all (-) eigenvalues","Maximum # of NCs 0 -> >0","Maximum # of NCs >0 -> >>0")
  col_names_combinedMetrics3 <- rbind("Minimum Total Abundance - Whole Community","Minimum Total Abundance - Probiotic Support Community","Minimum Total Abundance - Target Species","Minimum Total Abundance - Native Community","Minimum Relative Abundance - Probiotic Support Community","Minimum Relative Abundance - Target Species","Minimum Relative Abundance - Native Community","Minimum # of NC Target Species Diff > 0 (total abundance basis)","Minimum # of NC Target Species Diff < 0 (total abundance basis)","Minimum # of NC Target Species Diff > 0 (relative abundance basis)","Minimum # of NC Target Species Diff < 0 (relative abundance basis)","Minimum # of NC target species total abundance > 0.001","Minimum # of NC target species total abundance > 0.01","Minimum # of NC target species total abundance > 0.1","Minimum # of NC target species relative abundance > 0.001","Minimum # of NC target species relative abundance > 0.01","Minimum # of NC target species relative abundance > 0.1","Maxiumum Average Final Combined Community Evenness","Minimum # of NC with all (-) eigenvalues","Minimum # of NCs 0 -> >0","Minimum # of NCs >0 -> >>0")
  col_names_combinedMetrics4 <- rbind("Single Target Species Inoculation - # of NC target species total abundance > 0.001","Single Target Species Inoculation - # of NC target species total abundance > 0.01","Single Target Species Inoculation - # of NC target species total abundance > 0.1","Single Target Species Inoculation - # of NC target species relative abundance > 0.001","Single Target Species Inoculation - # of NC target species relative abundance > 0.01","Single Target Species Inoculation - # of NC target species relative abundance > 0.1","Single Target Speices Inoculation - Mean Final Community Evenness","Single Target Species Inoculation - Mean # of NC with all (-) eigenvalues","Single Target Species Inoculation - Mean target species total abundance","Single Target Species Inoculation - Mean Whole Community Total Abundance","Single Target Species Inoculation - Mean Native Community Total Abundance","Single Target Species Inoculation - Mean Target Species Relative Abundance","Single Target Species Inoculation - Mean Native Community Relative Abundance")
  col_names_combinedMetrics <- rbind(col_names_combinedMetrics1, col_names_combinedMetrics2, col_names_combinedMetrics3, col_names_combinedMetrics4)
  colnames(combinedMetrics) <- col_names_combinedMetrics
  combinedMetrics[1, 1] <- topPSC
  combinedMetrics[1, 2:20] <- averageMetrics[1:19]
  combinedMetrics[1, 21] <- averageMetrics[21]
  combinedMetrics[1, 22] <- averageMetrics[22]
  combinedMetrics[1, 23:41] <- maxMetrics[1:19]
  combinedMetrics[1, 42] <- maxMetrics[21]
  combinedMetrics[1, 43] <- maxMetrics[22]
  combinedMetrics[1, 44:62] <- minMetrics[1:19]
  combinedMetrics[1, 63] <- minMetrics[21]
  combinedMetrics[1, 64] <- minMetrics[22]
  combinedMetrics[1, 65:77] <- SingleKillerMetrics
  
  #writing KillerCommMetrics table into file_number directory
  csvFileName <- paste("KillerCommMetrics_file_number_",topPSC,".csv", sep="")
  csvFilePath <- paste0(wrkdir, "evolution_of_PSC_", PSC_number, "/Round_", z, "_data/",csvFileName) 
  write.table(KillerCommMetrics, file = csvFilePath, sep = ",", col.names = TRUE, row.names = F)
  
  #writing SingleKillerMetrics table into file_number directory
  csvFileName <- paste0("CombinedMetrics_file_number_",topPSC,".csv")
  csvFilePath <- paste0(wrkdir, "evolution_of_PSC_", PSC_number, "/Round_", z, "_data/",csvFileName) 
  write.table(combinedMetrics, file = csvFilePath, sep = ",", col.names = TRUE, row.names = F)
  
  TopKillerComm <- as.data.frame(KillerCommMetrics) %>% 
    arrange(desc(.[[22]])) %>%
    slice(1)
  
  topPSC <- as.numeric(TopKillerComm[, 20])
  
  csvFilePath <- paste0(wrkdir, "evolution_of_PSC_", PSC_number, "/Round_", z, "_data/KillerComm_stode_community_parameters/KillerComm_", topPSC, "_stode_final_parameters.csv")
  file.copy(csvFilePath, paste0(wrkdir, "evolution_of_PSC_", PSC_number, "/Round_", z, "_data"))
  
  csvFilePath <- paste0(wrkdir, "evolution_of_PSC_", PSC_number, "/Round_", z, "_data/KillerComm_stode_community_eigs/KillerComm_", topPSC, "_stode_final_eigs.csv")
  file.copy(csvFilePath, paste0(wrkdir, "evolution_of_PSC_", PSC_number, "/Round_", z, "_data"))
  
  csvFilePath <- paste0(wrkdir, "evolution_of_PSC_", PSC_number, "/Round_", z, "_data/KillerComm_stode_community_states/KillerComm_", topPSC, "_stode_final_states.csv")
  file.copy(csvFilePath, paste0(wrkdir, "evolution_of_PSC_", PSC_number, "/Round_", z, "_data"))
  
  print(paste0("finished round ", z))
  print(paste0("new best PSC = ", topPSC))
  print(paste0("number of parameters identified = ", nrow(newParms)))
  
  z <- z + 1
  if(all(parmcheck %in% oldcheck)) {
    break
  }
}
