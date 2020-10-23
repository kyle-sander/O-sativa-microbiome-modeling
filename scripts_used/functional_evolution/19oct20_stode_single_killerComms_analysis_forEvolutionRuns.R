  library(rootSolve)
  
  rm(list=ls())
  #Fill These Out
  num_of_comms <-1000
  file_number<-255
  num_of_KillerComm <- 2000
  num_of_species <- 11
  # mean and standard deviation for interaction coefficients  
  param_mean <- (-0.375)
  param_sd <- 1
  # mean and standard deviation for growth rates
  mu_mean <- 0.0002
  mu_sd <- 0.0002
  # assigns initial conditions, in this case all species start at same abundance
  start_amt<-0.001
  
  #initializing matricies and variables for cycling
  final_parameters <- matrix(ncol=30, nrow=num_of_comms) #might not be needed or used
  final_states <- matrix(ncol=5,nrow=num_of_comms) #might not be needed or used
  final_eigs <- matrix(ncol=5,nrow=num_of_comms) #might not be needed or used
  final_parameters_KillerComm <- matrix(ncol = (num_of_species*num_of_species+num_of_species), nrow = num_of_comms)
  final_states_KillerComm <- matrix(ncol = num_of_species, nrow = num_of_comms)
  final_eigs_KillerComm <- matrix(ncol=num_of_species, nrow=num_of_comms)
  final_parameters_singlekiller <- matrix(ncol = 42, nrow = num_of_comms)
  final_states_singlekiller <- matrix(ncol = 6, nrow = num_of_comms)
  final_eigs_singlekiller <- matrix(ncol=6,nrow=num_of_comms)
  n<-1 #begin loop counter 
  set.seed(file_number) #set seed for rnorm random sampling from normal distributions - setting seed with file_number variable, which should be unique for each set
  
  #reading in data from native community generation
  NativeCommunities_stode_final_parameters <- read.csv("~/Documents/SynCom_Modeling/NativeCommunities_stode2_final_parameters.csv", header=FALSE)
  NativeCommunities_stode_final_states <- read.csv("~/Documents/SynCom_Modeling/NativeCommunities_stode2_final_states.csv", header=FALSE)
  SKparamsForSKreassess<-read.csv("~/Documents/SynCom_Modeling/4/SingleKiller_stode_final_parameters.csv", header=FALSE) #reading in parameters from STSI inoculation, set 4 for import into STSI inoculation - reassessment to ensure abundance score = 430
  
  #begin single killer community generation
  for (n in 1:num_of_comms){
    muX1<-NativeCommunities_stode_final_parameters[n,1]
    muX2<-NativeCommunities_stode_final_parameters[n,2] 
    muX3<-NativeCommunities_stode_final_parameters[n,3]
    muX4<-NativeCommunities_stode_final_parameters[n,4]
    muX5<-NativeCommunities_stode_final_parameters[n,5]
    a11<-NativeCommunities_stode_final_parameters[n,6]
    a12<-NativeCommunities_stode_final_parameters[n,7]
    a13<-NativeCommunities_stode_final_parameters[n,8]
    a14<-NativeCommunities_stode_final_parameters[n,9]
    a15<-NativeCommunities_stode_final_parameters[n,10]
    a21<-NativeCommunities_stode_final_parameters[n,11]
    a22<-NativeCommunities_stode_final_parameters[n,12]
    a23<-NativeCommunities_stode_final_parameters[n,13]
    a24<-NativeCommunities_stode_final_parameters[n,14]
    a25<-NativeCommunities_stode_final_parameters[n,15]
    a31<-NativeCommunities_stode_final_parameters[n,16]
    a32<-NativeCommunities_stode_final_parameters[n,17]
    a33<-NativeCommunities_stode_final_parameters[n,18]
    a34<-NativeCommunities_stode_final_parameters[n,19]
    a35<-NativeCommunities_stode_final_parameters[n,20]
    a41<-NativeCommunities_stode_final_parameters[n,21]
    a42<-NativeCommunities_stode_final_parameters[n,22]
    a43<-NativeCommunities_stode_final_parameters[n,23]
    a44<-NativeCommunities_stode_final_parameters[n,24]
    a45<-NativeCommunities_stode_final_parameters[n,25]
    a51<-NativeCommunities_stode_final_parameters[n,26]
    a52<-NativeCommunities_stode_final_parameters[n,27]
    a53<-NativeCommunities_stode_final_parameters[n,28]
    a54<-NativeCommunities_stode_final_parameters[n,29]
    a55<-NativeCommunities_stode_final_parameters[n,30]
    muX6<-SKparamsForSKreassess[n,6]
    a16<-SKparamsForSKreassess[n,12]
    a61<-SKparamsForSKreassess[n,37]
    a26<-SKparamsForSKreassess[n,18]
    a62<-SKparamsForSKreassess[n,38]
    a36<-SKparamsForSKreassess[n,24]
    a63<-SKparamsForSKreassess[n,39]
    a46<-SKparamsForSKreassess[n,30]
    a64<-SKparamsForSKreassess[n,40]
    a56<-SKparamsForSKreassess[n,36]
    a65<-SKparamsForSKreassess[n,41]
    a66<-SKparamsForSKreassess[n,42]
    
    state<-c(X1=NativeCommunities_stode_final_states[n,1],X2=NativeCommunities_stode_final_states[n,2],X3=NativeCommunities_stode_final_states[n,3],X4=NativeCommunities_stode_final_states[n,4],X5=NativeCommunities_stode_final_states[n,5],X6=start_amt)
    parametersSingleKiller<-c(muX1,
                              muX2, 
                              muX3,
                              muX4,
                              muX5,
                              muX6,
                              a11,
                              a12,
                              a13,
                              a14,
                              a15,
                              a16,
                              a21,
                              a22,
                              a23,
                              a24,
                              a25,
                              a26,
                              a31,
                              a32,
                              a33,
                              a34,
                              a35,
                              a36,
                              a41,
                              a42,
                              a43,
                              a44,
                              a45,
                              a46,
                              a51,
                              a52,
                              a53,
                              a54,
                              a55,
                              a56,
                              a61,
                              a62,
                              a63,
                              a64,
                              a65,
                              a66)
    # function for generating differential equations according to the Lotka-Volterra model
    NonLinearGLV<-function(t, state, parameters) {
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
    out_4<-stode(y=state,fun=NonLinearGLV,parms=parametersSingleKiller, pos=TRUE)
    jac<-jacobian.full(y=c(out_4$y[1],out_4$y[2],out_4$y[3],out_4$y[4],out_4$y[5],out_4$y[6]),func=NonLinearGLV, parms=parametersSingleKiller)
    eig<-eigen(jac)
    final_states_singlekiller[n,]=out_4$y
    final_parameters_singlekiller[n,]<-parametersSingleKiller
    final_eigs_singlekiller[n,]<-Re(eig$values)
    n=n+1
  }
  dir.create(paste("~/Documents/SynCom_Modeling/",file_number,sep=""))
  dir.create(paste("~/Documents/SynCom_Modeling/",file_number,"/KillerComm_stode_community_eigs",sep=""))
  dir.create(paste("~/Documents/SynCom_Modeling/",file_number,"/KillerComm_stode_community_states",sep=""))
  dir.create(paste("~/Documents/SynCom_Modeling/",file_number,"/KillerComm_stode_community_parameters",sep=""))
  
  #copy native community files into the created directory
  #copying NC final parameters into set directory
  file.copy("~/Documents/SynCom_Modeling/NativeCommunities_stode_final_parameters.csv",paste("~/Documents/SynCom_Modeling/",file_number,sep=""))
  #copying NC final states into set directory
  file.copy("~/Documents/SynCom_Modeling/NativeCommunities_stode_final_states.csv",paste("~/Documents/SynCom_Modeling/",file_number,sep=""))
  #copying NC final eigs into set directory
  file.copy("~/Documents/SynCom_Modeling/NativeCommunities_stode_final_eigs.csv",paste("~/Documents/SynCom_Modeling/",file_number,sep=""))
  
  #writing single killer parameters file into directory named as file_number
  csvFileName <- paste("SingleKiller_stode_final_parameters.csv")
  csvFilePathSKparams <-paste("~/Documents/SynCom_Modeling/",file_number,"/",csvFileName,sep="") 
  write.table(final_parameters_singlekiller, file = csvFilePathSKparams, sep = ",", col.names = F, row.names = F)
  
  #writing single killer states file into directory named as file_number
  csvFileName <- paste("SingleKiller_stode_final_states.csv")
  csvFilePathSKstates <-paste("~/Documents/SynCom_Modeling/",file_number,"/",csvFileName,sep="") 
  write.table(final_states_singlekiller, file = csvFilePathSKstates, sep = ",", col.names = F, row.names = F)
  
  #writing single killer states file into directory named as file_number
  csvFileName <- paste("SingleKiller_stode_final_eigs.csv")
  csvFilePathSKeigs <-paste("~/Documents/SynCom_Modeling/",file_number,"/",csvFileName,sep="") 
  write.table(final_eigs_singlekiller, file = csvFilePathSKeigs, sep = ",", col.names = F, row.names = F)
  
  SingleKiller_stode_final_parameters <- read.csv(csvFilePathSKparams, header=FALSE)
  
  #begin KillerComm community generation, reset loop counters
  n <- 1  # loop counter for number of native + killer combos
  p <- 1  # loop counter for number of KillerComms
  
  # begins loop that will be repeated for each killercomm
  for (p in 1:num_of_KillerComm){
    #KillerComm growth rates
    muX7<-abs(rnorm(1, mean=mu_mean, sd=mu_sd))
    muX8<-abs(rnorm(1, mean=mu_mean, sd=mu_sd))
    muX9<-abs(rnorm(1, mean=mu_mean, sd=mu_sd))
    muX10<-abs(rnorm(1, mean=mu_mean, sd=mu_sd))
    muX11<-abs(rnorm(1, mean=mu_mean, sd=mu_sd))
    #killerComm interaction parameters with Killer
    a67<-rnorm(1,mean=param_mean,sd=param_sd)
    a68<-rnorm(1,mean=param_mean,sd=param_sd)
    a69<-rnorm(1,mean=param_mean,sd=param_sd)
    a610<-rnorm(1,mean=param_mean,sd=param_sd)
    a611<-rnorm(1,mean=param_mean,sd=param_sd)
    a116<-rnorm(1,mean=param_mean,sd=param_sd)
    a106<-rnorm(1,mean=param_mean,sd=param_sd)
    a96<-rnorm(1,mean=param_mean,sd=param_sd)
    a86<-rnorm(1,mean=param_mean,sd=param_sd)
    a76<-rnorm(1,mean=param_mean,sd=param_sd)
    # NativeComm Interactions with Killer Comm
    a17<-rnorm(1,mean=param_mean,sd=param_sd)
    a18<-rnorm(1,mean=param_mean,sd=param_sd)
    a19<-rnorm(1,mean=param_mean,sd=param_sd)
    a110<-rnorm(1,mean=param_mean,sd=param_sd)
    a1_11<-rnorm(1,mean=param_mean,sd=param_sd)
    a27<-rnorm(1,mean=param_mean,sd=param_sd)
    a28<-rnorm(1,mean=param_mean,sd=param_sd)
    a29<-rnorm(1,mean=param_mean,sd=param_sd)
    a210<-rnorm(1,mean=param_mean,sd=param_sd)
    a211<-rnorm(1,mean=param_mean,sd=param_sd)
    a37<-rnorm(1,mean=param_mean,sd=param_sd)
    a38<-rnorm(1,mean=param_mean,sd=param_sd)
    a39<-rnorm(1,mean=param_mean,sd=param_sd)
    a310<-rnorm(1,mean=param_mean,sd=param_sd)
    a311<-rnorm(1,mean=param_mean,sd=param_sd)
    a47<-rnorm(1,mean=param_mean,sd=param_sd)
    a48<-rnorm(1,mean=param_mean,sd=param_sd)
    a49<-rnorm(1,mean=param_mean,sd=param_sd)
    a410<-rnorm(1,mean=param_mean,sd=param_sd)
    a411<-rnorm(1,mean=param_mean,sd=param_sd)
    a57<-rnorm(1,mean=param_mean,sd=param_sd)
    a58<-rnorm(1,mean=param_mean,sd=param_sd)
    a59<-rnorm(1,mean=param_mean,sd=param_sd)
    a510<-rnorm(1,mean=param_mean,sd=param_sd)
    a511<-rnorm(1,mean=param_mean,sd=param_sd)
    #KillerComm interactions with native community
    a71<-rnorm(1,mean=param_mean,sd=param_sd)
    a72<-rnorm(1,mean=param_mean,sd=param_sd)
    a73<-rnorm(1,mean=param_mean,sd=param_sd)
    a74<-rnorm(1,mean=param_mean,sd=param_sd)
    a75<-rnorm(1,mean=param_mean,sd=param_sd)
    a77<-abs(rnorm(1,mean=param_mean,sd=param_sd))*-1
    a78<-rnorm(1,mean=param_mean,sd=param_sd)
    a79<-rnorm(1,mean=param_mean,sd=param_sd)
    a710<-rnorm(1,mean=param_mean,sd=param_sd)
    a711<-rnorm(1,mean=param_mean,sd=param_sd)
    a81<-rnorm(1,mean=param_mean,sd=param_sd)
    a82<-rnorm(1,mean=param_mean,sd=param_sd)
    a83<-rnorm(1,mean=param_mean,sd=param_sd)
    a84<-rnorm(1,mean=param_mean,sd=param_sd)
    a85<-rnorm(1,mean=param_mean,sd=param_sd)
    a87<-rnorm(1,mean=param_mean,sd=param_sd)
    a88<-abs(rnorm(1,mean=param_mean,sd=param_sd))*-1
    a89<-rnorm(1,mean=param_mean,sd=param_sd)
    a810<-rnorm(1,mean=param_mean,sd=param_sd)
    a811<-rnorm(1,mean=param_mean,sd=param_sd)
    a91<-rnorm(1,mean=param_mean,sd=param_sd)
    a92<-rnorm(1,mean=param_mean,sd=param_sd)
    a93<-rnorm(1,mean=param_mean,sd=param_sd)
    a94<-rnorm(1,mean=param_mean,sd=param_sd)
    a95<-rnorm(1,mean=param_mean,sd=param_sd)
    a97<-rnorm(1,mean=param_mean,sd=param_sd)
    a98<-rnorm(1,mean=param_mean,sd=param_sd)
    a99<-abs(rnorm(1,mean=param_mean,sd=param_sd))*-1
    a910<-rnorm(1,mean=param_mean,sd=param_sd)
    a911<-rnorm(1,mean=param_mean,sd=param_sd)
    a101<-rnorm(1,mean=param_mean,sd=param_sd)
    a102<-rnorm(1,mean=param_mean,sd=param_sd)
    a103<-rnorm(1,mean=param_mean,sd=param_sd)
    a104<-rnorm(1,mean=param_mean,sd=param_sd)
    a105<-rnorm(1,mean=param_mean,sd=param_sd)
    a107<-rnorm(1,mean=param_mean,sd=param_sd)
    a108<-rnorm(1,mean=param_mean,sd=param_sd)
    a109<-rnorm(1,mean=param_mean,sd=param_sd)
    a1010<-abs(rnorm(1,mean=param_mean,sd=param_sd))*-1
    a1011<-rnorm(1,mean=param_mean,sd=param_sd)
    a11_1<-rnorm(1,mean=param_mean,sd=param_sd)
    a112<-rnorm(1,mean=param_mean,sd=param_sd)
    a113<-rnorm(1,mean=param_mean,sd=param_sd)
    a114<-rnorm(1,mean=param_mean,sd=param_sd)
    a115<-rnorm(1,mean=param_mean,sd=param_sd)
    a117<-rnorm(1,mean=param_mean,sd=param_sd)
    a118<-rnorm(1,mean=param_mean,sd=param_sd)
    a119<-rnorm(1,mean=param_mean,sd=param_sd)
    a1110<-rnorm(1,mean=param_mean,sd=param_sd)
    a1111<-abs(rnorm(1,mean=param_mean,sd=param_sd))*-1
    
    # begins loop that will be repeated for each combo of killer + native community
    for (n in 1:num_of_comms){
      # reading in parameters from single killer model and native communities
      muX1<-NativeCommunities_stode_final_parameters[n,1]
      muX2<-NativeCommunities_stode_final_parameters[n,2] 
      muX3<-NativeCommunities_stode_final_parameters[n,3]
      muX4<-NativeCommunities_stode_final_parameters[n,4]
      muX5<-NativeCommunities_stode_final_parameters[n,5]
      a11<-NativeCommunities_stode_final_parameters[n,6]
      a12<-NativeCommunities_stode_final_parameters[n,7]
      a13<-NativeCommunities_stode_final_parameters[n,8]
      a14<-NativeCommunities_stode_final_parameters[n,9]
      a15<-NativeCommunities_stode_final_parameters[n,10]
      a21<-NativeCommunities_stode_final_parameters[n,11]
      a22<-NativeCommunities_stode_final_parameters[n,12]
      a23<-NativeCommunities_stode_final_parameters[n,13]
      a24<-NativeCommunities_stode_final_parameters[n,14]
      a25<-NativeCommunities_stode_final_parameters[n,15]
      a31<-NativeCommunities_stode_final_parameters[n,16]
      a32<-NativeCommunities_stode_final_parameters[n,17]
      a33<-NativeCommunities_stode_final_parameters[n,18]
      a34<-NativeCommunities_stode_final_parameters[n,19]
      a35<-NativeCommunities_stode_final_parameters[n,20]
      a41<-NativeCommunities_stode_final_parameters[n,21]
      a42<-NativeCommunities_stode_final_parameters[n,22]
      a43<-NativeCommunities_stode_final_parameters[n,23]
      a44<-NativeCommunities_stode_final_parameters[n,24]
      a45<-NativeCommunities_stode_final_parameters[n,25]
      a51<-NativeCommunities_stode_final_parameters[n,26]
      a52<-NativeCommunities_stode_final_parameters[n,27]
      a53<-NativeCommunities_stode_final_parameters[n,28]
      a54<-NativeCommunities_stode_final_parameters[n,29]
      a55<-NativeCommunities_stode_final_parameters[n,30]
      #single killer parameters 
      muX6<-SingleKiller_stode_final_parameters[n,6]
      a16<-SingleKiller_stode_final_parameters[n,12]
      a61<-SingleKiller_stode_final_parameters[n,37]
      a26<-SingleKiller_stode_final_parameters[n,18]
      a62<-SingleKiller_stode_final_parameters[n,38]
      a36<-SingleKiller_stode_final_parameters[n,24]
      a63<-SingleKiller_stode_final_parameters[n,39]
      a46<-SingleKiller_stode_final_parameters[n,30]
      a64<-SingleKiller_stode_final_parameters[n,40]
      a56<-SingleKiller_stode_final_parameters[n,36]
      a65<-SingleKiller_stode_final_parameters[n,41]
      a66<-SingleKiller_stode_final_parameters[n,42] 
      
      #At=matrix(c(a11,a12,a13,a14,a15,a16,a17,a18,a19,a110,a111,a21,a22,a23,a24,a25,a26,a27,a28,a29,a210,a211,a31,a32,a33,a34,a35,a36,a37,a38,a39,a310,a311,a41,a42,a43,a44,a45,a46,a47,a48,a49,a410,a411,a51,a52,a53,a54,a55,a56,a57,a58,a59,a510,a511,a61,a62,a63,a64,a65,a66,a67,a68,a69,a610,a611,a71,a72,a73,a74,a75,a76,a77,a78,a79,a710,a711,a81,a82,a83,a84,a85,a86,a87,a88,a89,a810,a811,a91,a92,a93,a94,a95,a96,a97,a98,a99,a910,a911,a101,a102,a103,a104,a105,a106,a107,a108,a109,a1010,a1011,a111,a112,a113,a114,a115,a116,a117,a118,a119,a1110,a1111), nrow=11, ncol=11)
      #A=t(At)
      #b=matrix(c(-muX1,-muX2,-muX3,-muX4,-muX5,-muX6,-muX7,-muX8,-muX9,-muX10,-muX11),nrow=11,ncol=1)
      #x<-solve(A,b)
      #if(x[1]>0&x[2]>0&x[3]>0&x[4]>0&x[5]>0){
      #zerocheck<-A%*%x-b
      
      # sets starting abundance to the final states obtained from the previous models and adds killercomm abundances
      state<-c(X1=NativeCommunities_stode_final_states[n,1],X2=NativeCommunities_stode_final_states[n,2],X3=NativeCommunities_stode_final_states[n,3],X4=NativeCommunities_stode_final_states[n,4],X5=NativeCommunities_stode_final_states[n,5],X6=start_amt,X7=start_amt,X8=start_amt,X9=start_amt,X10=start_amt,X11=start_amt)
      # function for generating differential equations according to the Lotka-Volterra model
      # assigns growth rates for all species and interaction coefficients as parameters for nonlinear solver
      parametersKillerComm<-c(muX1,
                              muX2, 
                              muX3,
                              muX4,
                              muX5,
                              muX6,
                              muX7,
                              muX8,
                              muX9,
                              muX10,
                              muX11,
                              a11,
                              a12,
                              a13,
                              a14,
                              a15,
                              a16,
                              a17,
                              a18,
                              a19,
                              a110,
                              a1_11,
                              a21,
                              a22,
                              a23,
                              a24,
                              a25,
                              a26,
                              a27,
                              a28,
                              a29,
                              a210,
                              a211,
                              a31,
                              a32,
                              a33,
                              a34,
                              a35,
                              a36,
                              a37,
                              a38,
                              a39,
                              a310,
                              a311,
                              a41,
                              a42,
                              a43,
                              a44,
                              a45,
                              a46,
                              a47,
                              a48,
                              a49,
                              a410,
                              a411,
                              a51,
                              a52,
                              a53,
                              a54,
                              a55,
                              a56,
                              a57,
                              a58,
                              a59,
                              a510,
                              a511,
                              a61,
                              a62,
                              a63,
                              a64,
                              a65,
                              a66,
                              a67,
                              a68,
                              a69,
                              a610,
                              a611,
                              a71,
                              a72,
                              a73,
                              a74,
                              a75,
                              a76,
                              a77,
                              a78,
                              a79,
                              a710,
                              a711,
                              a81,
                              a82,
                              a83,
                              a84,
                              a85,
                              a86,
                              a87,
                              a88,
                              a89,
                              a810,
                              a811,
                              a91,
                              a92,
                              a93,
                              a94,
                              a95,
                              a96,
                              a97,
                              a98,
                              a99,
                              a910,
                              a911,
                              a101,
                              a102,
                              a103,
                              a104,
                              a105,
                              a106,
                              a107,
                              a108,
                              a109,
                              a1010,
                              a1011,
                              a11_1,
                              a112,
                              a113,
                              a114,
                              a115,
                              a116,
                              a117,
                              a118,
                              a119,
                              a1110,
                              a1111)
      NonLinearGLV<-function(t, state, parameters) {
        with(as.list(c(state, parameters)),{
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
          list(c(dX1,dX2,dX3,dX4,dX5,dX6,dX7,dX8,dX9,dX10,dX11))
        })
      }
      out_4<-stode(y=state,fun=NonLinearGLV,parms=parametersKillerComm, pos=TRUE)
      jac<-jacobian.full(y=c(out_4$y[1],out_4$y[2],out_4$y[3],out_4$y[4],out_4$y[5],out_4$y[6],out_4$y[7],out_4$y[8],out_4$y[9],out_4$y[10],out_4$y[11]),func=NonLinearGLV, parms=parametersKillerComm)
      eig<-eigen(jac)
      final_states_KillerComm[n,]=out_4$y
      final_parameters_KillerComm[n,]<-parametersKillerComm
      final_eigs_KillerComm[n,]<-Re(eig$values)
      n=n+1
    }
    csvFileName <- paste("KillerComm_",p,"_stode_final_parameters.csv",sep="")
    csvFilePath <-paste("~/Documents/SynCom_Modeling/",file_number,"/KillerComm_stode_community_parameters/",csvFileName,sep="") 
    write.table(final_parameters_KillerComm, file = csvFilePath, sep = ",", col.names = F, row.names = F, append = TRUE)
    csvFileName2 <- paste("KillerComm_",p,"_stode_final_states.csv",sep="")
    csvFilePath2 <-paste("~/Documents/SynCom_Modeling/",file_number,"/KillerComm_stode_community_states/",csvFileName2,sep="") 
    write.table(final_states_KillerComm, file = csvFilePath2, sep = ",", col.names = F, row.names = F, append = TRUE)
    csvFileName3 <- paste("KillerComm_",p,"_stode_final_eigs.csv",sep="")
    csvFilePath3 <-paste("~/Documents/SynCom_Modeling/",file_number,"/KillerComm_stode_community_eigs/",csvFileName3,sep="") 
    write.table(final_eigs_KillerComm, file = csvFilePath3, sep = ",", col.names = F, row.names = F, append = TRUE)
    p = p+1
    print(p)
    n = 1
  }
  
  #end PSC generation, begin PSC analysis
  
  KillerCommMetrics<-matrix(nrow=num_of_KillerComm, ncol=22)
  col_names_KillerCommMetrics<-rbind("Average of Total Absolute Community Abundance across NCs","Average of Probiotic Support Community Absolute Community Abundance across NCs","Average of Target Species Absolute Community Abundance across NCs","Average of NC Absolute Community Abundance across NCs","Average Probiotic Support Community Relative Community Abundance","Average of Target Species Relative Community Abundance","Average of NC Combined Relative Community Abundance","# of NC Target Species Diff > 0 (total abundance basis)","# of NC Target Species Diff < 0 (total abundance basis)","# of NC Target Species Diff > 0 (relative abundance basis)","# of NC Target Species Diff < 0 (relative abundance basis)","# of NC Target Species Total Abundance > 0.001","# of NC Target Species Total Abundance > 0.01","# of NC Target Species Total Abundance > 0.1","# of NC Target Species Relative Abundance > 0.001","# of NC Target Species Relative Abundance > 0.01","# of NC Target Species Relative Abundance > 0.1","Average mean evenness of final communities across all NCs","# of NC with all (-) eigenvalues","Probiotic Support Community Number","# of NC whereby target species abundance = 0 in single target speices inoculation AND target species abundance > 0 in PSC inoculation","# of NC whereby target species abundance > 0 in single target speices inoculation AND target species abundance > in PSC inoculation than in single target species inoculation")
  colnames(KillerCommMetrics)<-col_names_KillerCommMetrics
  SingleKillerMetrics<-matrix(ncol=13, nrow=1)
  
  #importing Single Killer Final States
  SingleKillerFinStates<-read.csv(paste("~/Documents/SynCom_Modeling/",file_number,"/SingleKiller_stode_final_states.csv",sep=""), header=FALSE)
  SingleKillerEigs<-read.csv(paste("~/Documents/SynCom_Modeling/",file_number,"/SingleKiller_stode_final_eigs.csv",sep=""), header=FALSE)
  sumTotalAbund<-rowSums(SingleKillerFinStates, na.rm=TRUE) #total single killer community abundance
  sumNativeAbund<-rowSums(SingleKillerFinStates[,1:5], na.rm=TRUE) #native community abundance from 
  SingleKillerFinStates<-cbind(SingleKillerFinStates, sumTotalAbund, sumNativeAbund)
  SingleKillerRelStates<-matrix(ncol=18, nrow=num_of_comms)
  SingleKillerRelStates[,1]=SingleKillerFinStates[,1]/SingleKillerFinStates[,7]
  SingleKillerRelStates[,2]=SingleKillerFinStates[,2]/SingleKillerFinStates[,7]
  SingleKillerRelStates[,3]=SingleKillerFinStates[,3]/SingleKillerFinStates[,7]
  SingleKillerRelStates[,4]=SingleKillerFinStates[,4]/SingleKillerFinStates[,7]
  SingleKillerRelStates[,5]=SingleKillerFinStates[,5]/SingleKillerFinStates[,7]
  SingleKillerRelStates[,6]=SingleKillerFinStates[,6]/SingleKillerFinStates[,7]
  SingleKillerMetrics[,1]<-sum(SingleKillerFinStates[,6]>0.001, na.rm=TRUE) # of native comms target species total abundance > 0.001
  SingleKillerMetrics[,2]<-sum(SingleKillerFinStates[,6] >0.01, na.rm=TRUE) # of native comms target species total abundance > 0.01
  SingleKillerMetrics[,3]<-sum(SingleKillerFinStates[,6] >0.1, na.rm=TRUE) # of native comms target species total abundance > 0.1
  SingleKillerMetrics[,4]<-sum(SingleKillerRelStates[,6] >0.001, na.rm=TRUE) # of native comms target species total abundance > 0.001
  SingleKillerMetrics[,5]<-sum(SingleKillerRelStates[,6] >0.01, na.rm=TRUE) # of native comms target species total abundance > 0.01
  SingleKillerMetrics[,6]<-sum(SingleKillerRelStates[,6] >0.1, na.rm=TRUE) # of native comms target species total abundance > 0.1
  
  
  #calculating evenness for each Single Killer native comm
  n<-1
  for(n in 1:6){
    SingleKillerRelStates[,6+n]<-SingleKillerRelStates[,n]*(log(SingleKillerRelStates[,n]))
  } #calculating p*ln(p) for each species
  SingleKillerRelStates[,13]<-rowSums(SingleKillerRelStates[,7:12], na.rm=TRUE) #summing p*ln(p) values across all 11 species
  SingleKillerRelStates[,14]<-(SingleKillerRelStates[,13])*-1 #multiply sum by -1
  SingleKillerRelStates[,15]<-SingleKillerRelStates[,14]/log(6) # divide by log(6) - column 15 is evenness for each community
  SingleKillerRelStates[,16]<- SingleKillerRelStates[,6] #rel abundance target species single killer
  SingleKillerRelStates[,17]<- rowSums(SingleKillerRelStates[,1:5], na.rm=TRUE) #rel abundance native community single killer
  SingleKillerMetrics[,7]<-mean(SingleKillerRelStates[,15], na.rm=TRUE) #mean evenness of single killer final equilibrated communities
  SingleKillerMetrics[,8]<-length(which(rowSums(SingleKillerEigs<0)==6))  #number of single killer communities with all negative eigenvalues
  SingleKillerMetrics[,9]<- mean(SingleKillerFinStates[,6], na.rm=TRUE) #single killer target species total abundance
  SingleKillerMetrics[,10]<- mean(SingleKillerFinStates[,7], na.rm=TRUE) #single killer total community total abundance
  SingleKillerMetrics[,11]<- mean(SingleKillerFinStates[,8], na.rm=TRUE) #single killer native community total abundance
  SingleKillerMetrics[,12]<- mean(SingleKillerRelStates[,6], na.rm=TRUE) #single killer target species relative abuundance
  SingleKillerMetrics[,13]<- mean(SingleKillerRelStates[,17], na.rm=TRUE) #single killer native community relative abundance
  
  
  p<-1
  for(p in 1:num_of_KillerComm){
    #read in final states for KillerComm = p
    csvFileName <- paste("KillerComm_",p,"_stode_final_states.csv",sep="")
    csvFilePath <-paste("~/Documents/SynCom_Modeling/",file_number,"/KillerComm_stode_community_states/",csvFileName,sep="") 
    commStates<-read.table(file = csvFilePath, sep=",")
    #read in final eigenvalues for KillerComm = p
    csvFileName2 <- paste("KillerComm_",p,"_stode_final_eigs.csv",sep="")
    csvFilePath2 <-paste("~/Documents/SynCom_Modeling/",file_number,"/KillerComm_stode_community_eigs/",csvFileName2,sep="") 
    commEigs<-read.table(file = csvFilePath2, sep=",")
    
    
    #make zero vectors to tack onto as columns to commStates
    if(num_of_species!=11){
      num_of_columns_to_add<-11-num_of_species
      zeroColumns<-matrix(0L, ncol=num_of_columns_to_add, nrow=num_of_comms)
      commStates<-cbind(commStates, zeroColumns)
    }
    
    #Begin community metric calculations
    sumTotalAbund<-rowSums(commStates, na.rm=TRUE) #summing total absolute abundance of each final community
    sumKillerCommTotalAbund<-rowSums(commStates[,7:11], na.rm=TRUE)
    targetTotalAbund<-commStates[,6]
    nativeCommTotalAbund<-rowSums(commStates[,1:5], na.rm=TRUE)
    KillerCommTotalDiffs<-commStates[,6]-SingleKillerFinStates[,6]
    commStates<-cbind(commStates, sumTotalAbund, sumKillerCommTotalAbund, targetTotalAbund, nativeCommTotalAbund, KillerCommTotalDiffs)
    KillerCommMetrics[p,1]<-mean(commStates[,12]) #average of Total Absolute Community Abundance into KillerCommMetrics
    KillerCommMetrics[p,2]<-mean(commStates[,13]) #average of KillerComm Absolute Community Abundance into KillerCommMetrics
    KillerCommMetrics[p,3]<-mean(commStates[,14]) #average of Target Species Absolute Community Abundance into KillerCommMetrics
    KillerCommMetrics[p,4]<-mean(commStates[,15]) #average of NativeComm Absolute Community Abundance into KillerCommMetrics
    #calculating relative abundances of each community member
    commStatesRel<-matrix(ncol=25, nrow=num_of_comms)
    commStatesRel[,1]<-commStates[,1]/commStates[,12]
    commStatesRel[,2]<-commStates[,2]/commStates[,12]
    commStatesRel[,3]<-commStates[,3]/commStates[,12]
    commStatesRel[,4]<-commStates[,4]/commStates[,12]
    commStatesRel[,5]<-commStates[,5]/commStates[,12]
    commStatesRel[,6]<-commStates[,6]/commStates[,12]
    commStatesRel[,7]<-commStates[,7]/commStates[,12]
    commStatesRel[,8]<-commStates[,8]/commStates[,12]
    commStatesRel[,9]<-commStates[,9]/commStates[,12]
    commStatesRel[,10]<-commStates[,10]/commStates[,12]
    commStatesRel[,11]<-commStates[,11]/commStates[,12]
    sumKillerCommRelAbund<-rowSums(commStatesRel[,7:11])
    targetRelAbund<-commStatesRel[,6]
    nativeCommRelAbund<-rowSums(commStatesRel[,1:5])
    KillerCommRelDiffs<-commStatesRel[,6]-SingleKillerRelStates[,6]
    #calculating evenness for each final comm
    
    n<-1
    for(n in 1:num_of_species){
      commStatesRel[,11+n]<-commStatesRel[,n]*(log(commStatesRel[,n]))
    } #calculating p*ln(p) for each species
    commStatesRel[,23]<-rowSums(commStatesRel[,12:22], na.rm=TRUE) #summing p*ln(p) values across all species 
    commStatesRel[,24]<-(commStatesRel[,23])*-1 #multiply sum by -1
    commStatesRel[,25]<-commStatesRel[,24]/log(num_of_species) # divide by log(11) - this is final comm evenness 
    
    #entering metrics into KillerCommMetrics
    commStatesRel<-cbind(commStatesRel, sumKillerCommRelAbund, targetRelAbund, nativeCommRelAbund, KillerCommRelDiffs)
    commStatesRel[is.na(commStatesRel)]<-0
    KillerCommMetrics[p,5]<-mean(commStatesRel[,26], na.rm=TRUE) #average of KillerComm Combined Relative Community Abundance into KillerCommMetrics
    KillerCommMetrics[p,6]<-mean(commStatesRel[,27], na.rm=TRUE) #average of target species Relative Community Abundance into KillerCommMetrics
    KillerCommMetrics[p,7]<-mean(commStatesRel[,28], na.rm=TRUE) #average of NativeComm Combined Relative Community Abundance into KillerCommMetrics
    KillerCommMetrics[p,8]<-sum(commStates[,16] >0, na.rm=TRUE)# of Native Comms KillerComm target species Diff > 0 (total)
    KillerCommMetrics[p,9]<-sum(commStates[,16] <0, na.rm=TRUE)# of Native Comms KillerComm target species Diff < 0 (total)
    KillerCommMetrics[p,10]<-sum(commStatesRel[,29] >0, na.rm=TRUE) # of Native Comms KillerComm target species Diff > 0 (relative)
    KillerCommMetrics[p,11]<-sum(commStatesRel[,29] <0, na.rm=TRUE) # of Native Comms KillerComm target species Diff < 0 (relative)
    KillerCommMetrics[p,12]<-sum(commStates[,6] >0.001, na.rm=TRUE) # of native comms target species total abundance > 0.001
    KillerCommMetrics[p,13]<-sum(commStates[,6] >0.01, na.rm=TRUE) # of native comms target species total abundance > 0.01
    KillerCommMetrics[p,14]<-sum(commStates[,6] >0.1, na.rm=TRUE) # of native comms target species total abundance > 0.1
    KillerCommMetrics[p,15]<-sum(commStatesRel[,6] >0.001, na.rm=TRUE) # of comms target species relative abundance > 0.001
    KillerCommMetrics[p,16]<-sum(commStatesRel[,6] >0.01, na.rm=TRUE) # of comms target species relative abundance > 0.01
    KillerCommMetrics[p,17]<-sum(commStatesRel[,6] >0.1, na.rm=TRUE) # of comms target species relative abundance > 0.1
    KillerCommMetrics[p,18]<-mean(commStatesRel[,25], na.rm=TRUE) #mean evenness of final community
    num_of_neg_eigs<-rowSums(commEigs<0)
    KillerCommMetrics[p,19]<-length(which(num_of_neg_eigs==num_of_species)) # number of communities with all negative eigenvalues
    KillerCommMetrics[p,20]<-p #records probiotic support community number for row 
    
    # num of NC [0 in single species inoc] -> [>0 in probiotic support community inoc]
    #1. [Single Killer Fin States - species 6] AND [final States Killer Comm - species 6] into same matrix, two columns
    perform_metrics_table<-matrix(nrow=num_of_comms, ncol=5)
    perform_metrics_table[,1]<-commStates[,6]
    perform_metrics_table[,2]<-SingleKillerFinStates[,6]
    #2. column 3, give a 1 if single killer Fin state = 0
    perform_metrics_table[,3]<-ifelse(perform_metrics_table[,2]==0,1,0)
    #3. column 4, give a 1 if [killercomm-singleinoc]>0
    perform_metrics_table[,4]<-ifelse(perform_metrics_table[,1]>perform_metrics_table[,2],1,0)
    #4. column 5, add two columns
    perform_metrics_table[,5]<-perform_metrics_table[,3]+perform_metrics_table[,4]
    #5. NC of [0] -> [>0] is sums(column 5 = 2) - number of rows whose column #5 is = 2
    NCs_zero_to_pos<-sum(perform_metrics_table[,5]==2)
    #6. Put this # of NC inot KillerCommMetrics[p,21]
    KillerCommMetrics[p,21]<-NCs_zero_to_pos # of NC whereby target species abundance = 0 in single target speices inoculation AND target species abundance > 0 in PSC inoculation
    
    
    perform_metrics_table2<-matrix(nrow=num_of_comms, ncol=5)
    perform_metrics_table2[,1]<-commStates[,6] #final states from PSC
    perform_metrics_table2[,2]<-SingleKillerFinStates[,6] #final states from single target species inoculation
    # of NC [>0 in single species inoc] -> [>>0 in probiotic support community inoc] 
    #1. [Single Killer Fin States - species 6] AND [final States Killer Comm - species 6] into same matrix, two columns
    #2. column 3, give a 1 if single killer Fin state > 0
    perform_metrics_table2[,3]<-ifelse(perform_metrics_table2[,2]>0,1,0)
    #3. column 4, give a 1 if [killercomm-singleinoc]>0
    perform_metrics_table2[,4]<-ifelse(perform_metrics_table2[,1]>perform_metrics_table2[,2],1,0)
    #4. column 5, add two columns
    perform_metrics_table2[,5]<-perform_metrics_table2[,3]+perform_metrics_table2[,4]
    #5. NC of [0] -> [>0] is sums(column 5 = 2) - number of rows whose column #5 is = 2
    NCs_pos_to_greater<-sum(perform_metrics_table2[,5]==2)
    #6. Put this # of NC inot ParameterSensitivity[p,22]
    KillerCommMetrics[p,22]<-NCs_pos_to_greater # of NC whereby target species abundance > 0 in single target speices inoculation AND target species abundance > in PSC inoculation than in single target species inoculation
    
    rm(perform_metrics_table)
    rm(perform_metrics_table2)
    rm(commStatesRel)
    p=p+1
    print(p)
  }
  
  #writing KillerCommMetrics table into file_number directory
  csvFileName <- paste("KillerCommMetrics_file_number_",file_number,".csv", sep="")
  csvFilePath <-paste("~/Documents/SynCom_Modeling/",file_number,"/",csvFileName,sep="") 
  write.table(KillerCommMetrics, file = csvFilePath, sep = ",", col.names = TRUE, row.names = F)
  
  #writing SingleKillerMetrics table into file_number directory
  #csvFileName <- paste("SingleKillerMetrics.csv")
  #csvFilePath <-paste("~/Documents/SynCom_Modeling/",file_number,"/",csvFileName,sep="") 
  #write.table(SingleKillerMetrics, file = csvFilePath, sep = ",", col.names = F, row.names = F, append = TRUE)
  
  averageMetrics<-colMeans(KillerCommMetrics, na.rm=TRUE)
  maxMetrics<-apply(KillerCommMetrics,2,max, na.rm=TRUE) #Finds max of each column in KillerCommMetrics matrix
  minMetrics<-apply(KillerCommMetrics,2,min, na.rm=TRUE) #Finds min of each column in KillerCommMetrics matrix
  
  #combining metrics tables, generating 'combinedMetrics' and writing to file directory for set
  combinedMetrics<-matrix(nrow=1,ncol=77)
  col_names_combinedMetrics1<-rbind("file_set number","Average Total Abundance – Whole Community","Average Total Abundance – Probiotic Support Community","Average Total Abundance – Target Species","Average Total Abundance – Native Community","AverageRelative Abundance – Probiotic Support Community","Average Relative Abundance – Target Species","Average Relative Abundance – Native Community","Average # of NC Target Species Diff > 0 (total abundance, average of all probiotic support communities)","Average # of NC Target Species Diff < 0 (total abundance, average of all probiotic support communities)","Average # of NC Target Species Diff > 0 (relative abundance, average of all probiotic support communities)","Average # of NC Target Species Diff < 0 (relative abundance, average of all probiotic support communities)","Average # of NC target species total abundance > 0.001","Average # of NC target species total abundance > 0.01","Average # of NC target species total abundance > 0.1","Average # of NC target species relative abundance > 0.001","Average # of NC target species relative abundance > 0.01","Average # of NC target species relative abundance > 0.1","Average of Average Community Evenness","Average # of NC with all (-) eigenvalues","Average # of NCs 0 -> >0","Average # of NCs >0 -> >>0")
  col_names_combinedMetrics2<-rbind("Maximum Total Abundance - Whole Community","Maximum Total Abundance - Probiotic Support Community","Maximum Total Abundance - Target Species","Maximum Total Abundance - Native Community","Maximum Relative Abundance - Probiotic Support Community","Maximum Relative Abundance - Target Species","Maximum Relative Abundance - Native Community","Maximum # of NC Target Species Diff > 0 (total abundance basis)","Maximum # of NC Target Species Diff < 0 (total abundance basis)","Maximum # of NC Target Species Diff > 0 (relative abundance basis)","Maximum # of NC Target Species Diff < 0 (relative abundance basis)","Maximum # of NC target species total abundance > 0.001","Maximum # of NC target species total abundance > 0.01","Maximum # of NC target species total abundance > 0.1","Maximum # of NC target species relative abundance > 0.001","Maximum # of NC target species relative abundance > 0.01","Maximum # of NC target species relative abundance > 0.1","Maxiumum Average Final Combined Community Evenness","Maximum # of NC with all (-) eigenvalues","Maximum # of NCs 0 -> >0","Maximum # of NCs >0 -> >>0")
  col_names_combinedMetrics3<-rbind("Minimum Total Abundance - Whole Community","Minimum Total Abundance - Probiotic Support Community","Minimum Total Abundance - Target Species","Minimum Total Abundance - Native Community","Minimum Relative Abundance - Probiotic Support Community","Minimum Relative Abundance - Target Species","Minimum Relative Abundance - Native Community","Minimum # of NC Target Species Diff > 0 (total abundance basis)","Minimum # of NC Target Species Diff < 0 (total abundance basis)","Minimum # of NC Target Species Diff > 0 (relative abundance basis)","Minimum # of NC Target Species Diff < 0 (relative abundance basis)","Minimum # of NC target species total abundance > 0.001","Minimum # of NC target species total abundance > 0.01","Minimum # of NC target species total abundance > 0.1","Minimum # of NC target species relative abundance > 0.001","Minimum # of NC target species relative abundance > 0.01","Minimum # of NC target species relative abundance > 0.1","Maxiumum Average Final Combined Community Evenness","Minimum # of NC with all (-) eigenvalues","Minimum # of NCs 0 -> >0","Minimum # of NCs >0 -> >>0")
  col_names_combinedMetrics4<-rbind("Single Target Species Inoculation - # of NC target species total abundance > 0.001","Single Target Species Inoculation - # of NC target species total abundance > 0.01","Single Target Species Inoculation - # of NC target species total abundance > 0.1","Single Target Species Inoculation - # of NC target species relative abundance > 0.001","Single Target Species Inoculation - # of NC target species relative abundance > 0.01","Single Target Species Inoculation - # of NC target species relative abundance > 0.1","Single Target Speices Inoculation - Mean Final Community Evenness","Single Target Species Inoculation - Mean # of NC with all (-) eigenvalues","Single Target Species Inoculation - Mean target species total abundance","Single Target Species Inoculation - Mean Whole Community Total Abundance","Single Target Species Inoculation - Mean Native Community Total Abundance","Single Target Species Inoculation - Mean Target Species Relative Abundance","Single Target Species Inoculation - Mean Native Community Relative Abundance")
  col_names_combinedMetrics<-rbind(col_names_combinedMetrics1,col_names_combinedMetrics2,col_names_combinedMetrics3,col_names_combinedMetrics4)
  colnames(combinedMetrics)<-col_names_combinedMetrics
  combinedMetrics[1,1]<-file_number
  combinedMetrics[1,2:20]<-averageMetrics[1:19]
  combinedMetrics[1,21]<-averageMetrics[21]
  combinedMetrics[1,22]<-averageMetrics[22]
  combinedMetrics[1,23:41]<-maxMetrics[1:19]
  combinedMetrics[1,42]<-maxMetrics[21]
  combinedMetrics[1,43]<-maxMetrics[22]
  combinedMetrics[1,44:62]<-minMetrics[1:19]
  combinedMetrics[1,63]<-minMetrics[21]
  combinedMetrics[1,64]<-minMetrics[22]
  combinedMetrics[1,65:77]<-SingleKillerMetrics
  
  #writing SingleKillerMetrics table into file_number directory
  csvFileName <- paste("CombinedMetrics_file_number_",file_number,".csv",sep="")
  csvFilePath <-paste("~/Documents/SynCom_Modeling/",file_number,"/",csvFileName,sep="") 
  write.table(combinedMetrics, file = csvFilePath, sep = ",", col.names = TRUE, row.names = F)
