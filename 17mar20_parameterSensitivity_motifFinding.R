      #install.packages("rootSolve")
      library(rootSolve)
      
      #this script imports data for high-performing communities parameter set(s) which were generated using community generating scripts 
      #1. Identify High Performing PSC's from KillerCommMetrics Output
      #2. ......
      #3. ......
      
      rm(list=ls())
      #Fill These Out
      num_of_comms <-1000
      num_of_species <- 11
      # assigns initial conditions, in thi  s case all species start at same abundance
      start_amt<-0.001
      file_number<-3 #file number name where KillerComm and Single Killer files are stored
      PSC_number<-957 #Probiotic Support Community identifier number - Probiotic support community to analyze common parameters
      num_of_Common_KillerCommParams<-80 #number of parameters to cycle through and, separately, set to zero and assess performance
      #stode_or_desolve<- #if using stode generated community, input 'stode', if using desolve generated communities, input 'desolve'
      
      #initializing matricies and variables for cycling - SOME OF THESE AREN'T NEEDED/USED
      final_parameters <- matrix(ncol=30, nrow=num_of_comms)
      final_states <- matrix(ncol=5,nrow=num_of_comms)
      final_eigs <- matrix(ncol=5,nrow=num_of_comms)
      final_parameters_KillerComm <- matrix(ncol = (num_of_species*num_of_species+num_of_species), nrow = num_of_comms)
      final_states_KillerComm <- matrix(ncol = num_of_species, nrow = num_of_comms)
      final_eigs_KillerComm <- matrix(ncol=num_of_species,nrow=num_of_comms)
      final_parameters_singlekiller <- matrix(ncol = 42, nrow = num_of_comms)
      final_states_singlekiller <- matrix(ncol = 6, nrow = num_of_comms)
      final_eigs_singlekiller <- matrix(ncol=6,nrow=num_of_comms)
      CommonKillerCommParams<-matrix(nrow=1, ncol=num_of_Common_KillerCommParams+10)
      ParameterSensitivity<-matrix(nrow=num_of_Common_KillerCommParams+1, ncol=11)
      
      #Importing in needed parameter and state data from simulations
      #reading in parameters and states from native community generation
      #reading in native community final parameters
      csvFilePath<-paste("~/Documents/SynCom_Modeling/Native_Communities/stodeNCs/NativeCommunities_stode_final_parameters.csv")
      NativeCommunities_stode_final_parameters<-read.table(file = csvFilePath, sep=",",header = FALSE)
      
      #NEED TO IMPORT NATIVE COMMUNITIES FINAL STATES
      csvFilePath<-paste("~/Documents/SynCom_Modeling/Native_Communities/stodeNCs/NativeCommunities_stode_final_states.csv")
      NativeCommunities_stode_final_states<-read.table(file = csvFilePath, sep=",",header = FALSE)
      
      #reading in selected PSC parameters 
      csvFilePath2<-paste("~/Documents/SynCom_Modeling/",file_number,"/KillerComm_stode_community_parameters/KillerComm_",PSC_number,"_stode_final_parameters.csv", sep="")
      KillerComm_stode_highPerformingParams <- read.csv(file = csvFilePath2, sep=",",header = FALSE)
      
      #reading in selected PSC states
      csvFilePath3<-paste("~/Documents/SynCom_Modeling/",file_number,"/KillerComm_stode_community_states/KillerComm_",PSC_number,"_stode_final_states.csv", sep="")
      KillerComm_stode_highPerformingStates<-read.csv(file = csvFilePath3, sep=",",header = FALSE)
      
      #reading in selected PSC eigs
      csvFilePath4<-paste("~/Documents/SynCom_Modeling/",file_number,"/KillerComm_stode_community_eigs/KillerComm_",PSC_number,"_stode_final_eigs.csv", sep="")
      KillerComm_stode_highPerformingEigs<-read.csv(file = csvFilePath4, sep=",",header = FALSE)
      
      #reading in single target species inoculation final states
      csvFilePath5<-paste("~/Documents/SynCom_Modeling/",file_number,"/SingleKiller_stode_final_states.csv", sep="")
      SingleKillerFinStates<-read.csv(file = csvFilePath5, sep=",",header = FALSE)
      
      #reading in single target species inoculation final parameters
      csvFilePath5<-paste("~/Documents/SynCom_Modeling/",file_number,"/SingleKiller_stode_final_parameters.csv", sep="")
      SingleKiller_stode_final_parameters<-read.csv(file = csvFilePath5, sep=",",header = FALSE)
      
      
      
      #Import KillerComm_stode_final_states as well to calculate baseline, no parameter change and insert into first row of table
      ParameterSensitivity[1,1]<-paste("baseline - no parameters omitted") #This is identifier for scenario when no parameters have been augmented/omitted, the 'baseline' - 
      #use this matrix position as row label to indicate as such - instead of "0"
      ParameterSensitivity[1,3]<-sum(KillerComm_stode_highPerformingStates[,6]>0.001, na.rm=TRUE) #this is calculating which communities are >0.001 right now
      #want to shift this to be diffs > 0 as needed parameter - its cool this is next column, need column headers 
      KillerCommTotalDiffs<-KillerComm_stode_highPerformingStates[,6]-SingleKillerFinStates[,6]
      ParameterSensitivity[1,4]<-sum(KillerCommTotalDiffs>0, na.rm=TRUE)# of Native Comms KillerComm target species Diff > 0 (total)
      
      #baseline [0]->[>0] and [>0]->[>>0] - read in from KillerCommMetrics......
      #writing KillerCommMetrics table into file_number directory
      csvFileName <- paste("KillerCommMetrics_file_number_",file_number,".csv", sep="")
      csvFilePath <-paste("~/Documents/SynCom_Modeling/",file_number,"/",csvFileName,sep="") 
      KillerCommMetrics<-read.table(file = csvFilePath, sep = ",", stringsAsFactors = F)
      
      ParameterSensitivity[1,8]<-as.numeric(KillerCommMetrics[PSC_number+1,21]) #reading in [0]->[>0] for baseline, from KillerCommMetrics
      ParameterSensitivity[1,10]<-as.numeric(KillerCommMetrics[PSC_number+1,22]) #reading in [>0]->[>>0] for baseline, from KillerCommMetrics
      
      col_names<-rbind("Parameter Number being omitted/set to 0","Original Parameter Value","# of NC's target total abundance > 0.001","# of NC target species Diff > 0 (total abundance)","#num of NC with all negative eigenvalues","% performance when parameter omitted","Parameter","# of NC target species [0-single target inoc] -> [>0 - PSC inoc]","% Performance Difference [0-single target inoc] -> [>0 - PSC inoc]","# of NC target species [>0 - single target inoc] -> [>>0 - PSC inoc] ","% Performance Difference [>0 - single target inoc] -> [>>0 - PSC inoc]")
      colnames(ParameterSensitivity)<-col_names
      
      num_of_neg_eigs<-rowSums(KillerComm_stode_highPerformingEigs<0)  
      ParameterSensitivity[1,5]<-length(which(num_of_neg_eigs==num_of_species))# num of (-) eigenvalues for single target species inoculation
      
      
      
      #Columns 6 and 7 = # of NC(zero -> something) AND # of NC(something -> something higher) - splitting the Diffs > 0 into two categories
      #NEED TO CALCULATE AND ADD THIS IN
      
      #begin KillerComm community generation, reset loop counters
      n <- 1  # loop counter for number of native + killer combos
      p <- 1  # loop counter for number of KillerComm Common Parameters individually set to zero and reassess
      
      # begins loop that will be repeated, setting each individual parameter of CommonKillerCommParams, separately, 
      #to zero and then re-assessing support community performance against all 1000 native communities
      for (p in 1:num_of_Common_KillerCommParams){
        #import parameters from KillerComms that are common to all Native Comms
        CommonKillerCommParams[1,1]<-KillerComm_stode_highPerformingParams[1000,18] #a17
        CommonKillerCommParams[1,2]<-KillerComm_stode_highPerformingParams[1000,19] #a18
        CommonKillerCommParams[1,3]<-KillerComm_stode_highPerformingParams[1000,20] #a19
        CommonKillerCommParams[1,4]<-KillerComm_stode_highPerformingParams[1000,21] #a110
        CommonKillerCommParams[1,5]<-KillerComm_stode_highPerformingParams[1000,22] #a1_11
        CommonKillerCommParams[1,6]<-KillerComm_stode_highPerformingParams[1000,29] #a27
        CommonKillerCommParams[1,7]<-KillerComm_stode_highPerformingParams[1000,30] #a28
        CommonKillerCommParams[1,8]<-KillerComm_stode_highPerformingParams[1000,31] #a29
        CommonKillerCommParams[1,9]<-KillerComm_stode_highPerformingParams[1000,32] #a210
        CommonKillerCommParams[1,10]<-KillerComm_stode_highPerformingParams[1000,33] #a211
        CommonKillerCommParams[1,11]<-KillerComm_stode_highPerformingParams[1000,40] #a37
        CommonKillerCommParams[1,12]<-KillerComm_stode_highPerformingParams[1000,41] #a38
        CommonKillerCommParams[1,13]<-KillerComm_stode_highPerformingParams[1000,42] #a39
        CommonKillerCommParams[1,14]<-KillerComm_stode_highPerformingParams[1000,43] #a310
        CommonKillerCommParams[1,15]<-KillerComm_stode_highPerformingParams[1000,44] #a311
        CommonKillerCommParams[1,16]<-KillerComm_stode_highPerformingParams[1000,51] #a47
        CommonKillerCommParams[1,17]<-KillerComm_stode_highPerformingParams[1000,52] #a48
        CommonKillerCommParams[1,18]<-KillerComm_stode_highPerformingParams[1000,53] #a49
        CommonKillerCommParams[1,19]<-KillerComm_stode_highPerformingParams[1000,54] #a410
        CommonKillerCommParams[1,20]<-KillerComm_stode_highPerformingParams[1000,55] #a411
        CommonKillerCommParams[1,21]<-KillerComm_stode_highPerformingParams[1000,62] #a57
        CommonKillerCommParams[1,22]<-KillerComm_stode_highPerformingParams[1000,63] #a58
        CommonKillerCommParams[1,23]<-KillerComm_stode_highPerformingParams[1000,64] #a59
        CommonKillerCommParams[1,24]<-KillerComm_stode_highPerformingParams[1000,65] #a510
        CommonKillerCommParams[1,25]<-KillerComm_stode_highPerformingParams[1000,66] #a511
        CommonKillerCommParams[1,26]<-KillerComm_stode_highPerformingParams[1000,73] #a67
        CommonKillerCommParams[1,27]<-KillerComm_stode_highPerformingParams[1000,74] #a68
        CommonKillerCommParams[1,28]<-KillerComm_stode_highPerformingParams[1000,75] #a69
        CommonKillerCommParams[1,29]<-KillerComm_stode_highPerformingParams[1000,76] #a610
        CommonKillerCommParams[1,30]<-KillerComm_stode_highPerformingParams[1000,77] #a611
        CommonKillerCommParams[1,31]<-KillerComm_stode_highPerformingParams[1000,78] #a71
        CommonKillerCommParams[1,32]<-KillerComm_stode_highPerformingParams[1000,79] #a72
        CommonKillerCommParams[1,33]<-KillerComm_stode_highPerformingParams[1000,80] #a73
        CommonKillerCommParams[1,34]<-KillerComm_stode_highPerformingParams[1000,81] #a74
        CommonKillerCommParams[1,35]<-KillerComm_stode_highPerformingParams[1000,82] #a75
        CommonKillerCommParams[1,36]<-KillerComm_stode_highPerformingParams[1000,83] #a76
        CommonKillerCommParams[1,37]<-KillerComm_stode_highPerformingParams[1000,85] #a78
        CommonKillerCommParams[1,38]<-KillerComm_stode_highPerformingParams[1000,86] #a79
        CommonKillerCommParams[1,39]<-KillerComm_stode_highPerformingParams[1000,87] #a710
        CommonKillerCommParams[1,40]<-KillerComm_stode_highPerformingParams[1000,88] #a711
        CommonKillerCommParams[1,41]<-KillerComm_stode_highPerformingParams[1000,89] #a81
        CommonKillerCommParams[1,42]<-KillerComm_stode_highPerformingParams[1000,90] #a82
        CommonKillerCommParams[1,43]<-KillerComm_stode_highPerformingParams[1000,91] #a83
        CommonKillerCommParams[1,44]<-KillerComm_stode_highPerformingParams[1000,92] #a84
        CommonKillerCommParams[1,45]<-KillerComm_stode_highPerformingParams[1000,93] #a85
        CommonKillerCommParams[1,46]<-KillerComm_stode_highPerformingParams[1000,94] #a86
        CommonKillerCommParams[1,47]<-KillerComm_stode_highPerformingParams[1000,95] #a87
        CommonKillerCommParams[1,48]<-KillerComm_stode_highPerformingParams[1000,97] #a89
        CommonKillerCommParams[1,49]<-KillerComm_stode_highPerformingParams[1000,98] #a810
        CommonKillerCommParams[1,50]<-KillerComm_stode_highPerformingParams[1000,99] #a811
        CommonKillerCommParams[1,51]<-KillerComm_stode_highPerformingParams[1000,100] #a91
        CommonKillerCommParams[1,52]<-KillerComm_stode_highPerformingParams[1000,101] #a92
        CommonKillerCommParams[1,53]<-KillerComm_stode_highPerformingParams[1000,102] #a93
        CommonKillerCommParams[1,54]<-KillerComm_stode_highPerformingParams[1000,103] #a94
        CommonKillerCommParams[1,55]<-KillerComm_stode_highPerformingParams[1000,104] #a95
        CommonKillerCommParams[1,56]<-KillerComm_stode_highPerformingParams[1000,105] #a96
        CommonKillerCommParams[1,57]<-KillerComm_stode_highPerformingParams[1000,106] #a97
        CommonKillerCommParams[1,58]<-KillerComm_stode_highPerformingParams[1000,107] #a98
        CommonKillerCommParams[1,59]<-KillerComm_stode_highPerformingParams[1000,109] #a910
        CommonKillerCommParams[1,60]<-KillerComm_stode_highPerformingParams[1000,110] #a911
        CommonKillerCommParams[1,61]<-KillerComm_stode_highPerformingParams[1000,111] #a101
        CommonKillerCommParams[1,62]<-KillerComm_stode_highPerformingParams[1000,112] #a102
        CommonKillerCommParams[1,63]<-KillerComm_stode_highPerformingParams[1000,113] #a103
        CommonKillerCommParams[1,64]<-KillerComm_stode_highPerformingParams[1000,114] #a104
        CommonKillerCommParams[1,65]<-KillerComm_stode_highPerformingParams[1000,115] #a105
        CommonKillerCommParams[1,66]<-KillerComm_stode_highPerformingParams[1000,116] #a106
        CommonKillerCommParams[1,67]<-KillerComm_stode_highPerformingParams[1000,117] #a107
        CommonKillerCommParams[1,68]<-KillerComm_stode_highPerformingParams[1000,118] #a108
        CommonKillerCommParams[1,69]<-KillerComm_stode_highPerformingParams[1000,119] #a109
        CommonKillerCommParams[1,70]<-KillerComm_stode_highPerformingParams[1000,121] #a1011
        CommonKillerCommParams[1,71]<-KillerComm_stode_highPerformingParams[1000,122] #a11_1
        CommonKillerCommParams[1,72]<-KillerComm_stode_highPerformingParams[1000,123] #a112
        CommonKillerCommParams[1,73]<-KillerComm_stode_highPerformingParams[1000,124] #a113
        CommonKillerCommParams[1,74]<-KillerComm_stode_highPerformingParams[1000,125] #a114
        CommonKillerCommParams[1,75]<-KillerComm_stode_highPerformingParams[1000,126] #a115
        CommonKillerCommParams[1,76]<-KillerComm_stode_highPerformingParams[1000,127] #a116
        CommonKillerCommParams[1,77]<-KillerComm_stode_highPerformingParams[1000,128] #a117
        CommonKillerCommParams[1,78]<-KillerComm_stode_highPerformingParams[1000,129] #a118
        CommonKillerCommParams[1,79]<-KillerComm_stode_highPerformingParams[1000,130] #a119
        CommonKillerCommParams[1,80]<-KillerComm_stode_highPerformingParams[1000,131] #a1110
       
       #set parameter number p from CommonKillerCommParams equal to zero
        ParameterSensitivity[p+1,1]<-p #denotes which index of p is currently running
        ParameterSensitivity[p+1,2]<-CommonKillerCommParams[1,p] #denotes original 
        CommonKillerCommParams[1,p]<-0 #remember to set this back to original number prior to - fixed by reinitializing CommonKillerCommParams each cycle of p
        #this isnt going to work, because when p = 73 it is going to set nothing equal to 0 and then 
        #import Common KillerComm Params 
        
        #KillerComm growth rates
        muX7<-KillerComm_stode_highPerformingParams[1000,7] #muX7
        muX8<-KillerComm_stode_highPerformingParams[1000,8] #muX8
        muX9<-KillerComm_stode_highPerformingParams[1000,9] #muX9
        muX10<-KillerComm_stode_highPerformingParams[1000,10] #muX10
        muX11<-KillerComm_stode_highPerformingParams[1000,11] #muX11
        #killerComm interaction parameters with Killer
        a67<-CommonKillerCommParams[1,26] #a67
        a68<-CommonKillerCommParams[1,27] #a68
        a69<-CommonKillerCommParams[1,28] #a69
        a610<-CommonKillerCommParams[1,29] #a610
        a611<-CommonKillerCommParams[1,30] #a611
        a116<-CommonKillerCommParams[1,76] #a116
        a106<-CommonKillerCommParams[1,66] #a106
        a96<-CommonKillerCommParams[1,56] #a96
        a86<-CommonKillerCommParams[1,46] #a86
        a76<-CommonKillerCommParams[1,36] #a76
        # NativeComm Interactions with Killer Comm
        a17<-CommonKillerCommParams[1,1] #a17
        a18<-CommonKillerCommParams[1,2] #a18
        a19<-CommonKillerCommParams[1,3] #a19
        a110<-CommonKillerCommParams[1,4] #a110
        a1_11<-CommonKillerCommParams[1,5] #a1_11
        a27<-CommonKillerCommParams[1,6] #a27
        a28<-CommonKillerCommParams[1,7] #a28
        a29<-CommonKillerCommParams[1,8] #a29
        a210<-CommonKillerCommParams[1,9] #a210
        a211<-CommonKillerCommParams[1,10] #a211
        a37<-CommonKillerCommParams[1,11] #a37
        a38<-CommonKillerCommParams[1,12] #a38
        a39<-CommonKillerCommParams[1,13] #a39
        a310<-CommonKillerCommParams[1,14] #a310
        a311<-CommonKillerCommParams[1,15] #a311
        a47<-CommonKillerCommParams[1,16] #a47
        a48<-CommonKillerCommParams[1,17] #a48
        a49<-CommonKillerCommParams[1,18] #a49
        a410<-CommonKillerCommParams[1,19] #a410
        a411<-CommonKillerCommParams[1,20] #a411
        a57<-CommonKillerCommParams[1,21] #a57
        a58<-CommonKillerCommParams[1,22] #a58
        a59<-CommonKillerCommParams[1,23] #a59
        a510<-CommonKillerCommParams[1,24] #a510
        a511<-CommonKillerCommParams[1,25] #a511
        #KillerComm interactions with native community
        a71<-CommonKillerCommParams[1,31] #a71
        a72<-CommonKillerCommParams[1,32] #a72
        a73<-CommonKillerCommParams[1,33] #a73
        a74<-CommonKillerCommParams[1,34] #a74
        a75<-CommonKillerCommParams[1,35] #a75
        a77<-KillerComm_stode_highPerformingParams[1000,84] #a77
        a78<-CommonKillerCommParams[1,37] #a78
        a79<-CommonKillerCommParams[1,38] #a79
        a710<-CommonKillerCommParams[1,39] #a710
        a711<-CommonKillerCommParams[1,40] #a711
        a81<-CommonKillerCommParams[1,41] #a81
        a82<-CommonKillerCommParams[1,42] #a82
        a83<-CommonKillerCommParams[1,43] #a83
        a84<-CommonKillerCommParams[1,44] #a84
        a85<-CommonKillerCommParams[1,45] #a85
        a87<-CommonKillerCommParams[1,47] #a87
        a88<-KillerComm_stode_highPerformingParams[1000,96] #a88
        a89<-CommonKillerCommParams[1,48] #a89
        a810<-CommonKillerCommParams[1,49] #a810
        a811<-CommonKillerCommParams[1,50] #a811
        a91<-CommonKillerCommParams[1,51] #a91
        a92<-CommonKillerCommParams[1,52] #a92
        a93<-CommonKillerCommParams[1,53] #a93
        a94<-CommonKillerCommParams[1,54] #a94
        a95<-CommonKillerCommParams[1,55] #a95
        a97<-CommonKillerCommParams[1,57] #a97
        a98<-CommonKillerCommParams[1,58] #a98
        a99<-KillerComm_stode_highPerformingParams[1000,108] #a99
        a910<-CommonKillerCommParams[1,59] #a910
        a911<-CommonKillerCommParams[1,60] #a911
        a101<-CommonKillerCommParams[1,61] #a101
        a102<-CommonKillerCommParams[1,62] #a102
        a103<-CommonKillerCommParams[1,63] #a103
        a104<-CommonKillerCommParams[1,64] #a104
        a105<-CommonKillerCommParams[1,65] #a105
        a107<-CommonKillerCommParams[1,67] #a107
        a108<-CommonKillerCommParams[1,68] #a108
        a109<-CommonKillerCommParams[1,69] #a109
        a1010<-KillerComm_stode_highPerformingParams[1000,120] #a1010
        a1011<-CommonKillerCommParams[1,70] #a1011
        a11_1<-CommonKillerCommParams[1,71] #a11_1
        a112<-CommonKillerCommParams[1,72] #a112
        a113<-CommonKillerCommParams[1,73] #a113
        a114<-CommonKillerCommParams[1,74] #a114
        a115<-CommonKillerCommParams[1,75] #a115
        a117<-CommonKillerCommParams[1,77] #a117
        a118<-CommonKillerCommParams[1,78] #a118
        a119<-CommonKillerCommParams[1,79] #a119
        a1110<-CommonKillerCommParams[1,80] #a1110
        a1111<-KillerComm_stode_highPerformingParams[1000,132] #a1111
        
        # begins loop that will be repeated for each combo of killer + native community
        for (n in 1:num_of_comms){
          #import 
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
      
        #make column 2 the original parameter value
        ParameterSensitivity[p+1,3]<-sum(final_states_KillerComm[,6] > 0.001, na.rm=TRUE) #num of native comms target total abundance > 0.001
        
        # num of native comms target total abundance diff from single killer > 0
        KillerCommTotalDiffs<-final_states_KillerComm[,6]-SingleKillerFinStates[,6]
        ParameterSensitivity[p+1,4]<-sum(KillerCommTotalDiffs>0, na.rm=TRUE)# of Native Comms KillerComm target species Diff > 0 (total)
        
        # num of native comms with all Re(eigs) < 0  
        num_of_neg_eigs<-rowSums(final_eigs_KillerComm<0)  
        ParameterSensitivity[p+1,5]<-length(which(num_of_neg_eigs==num_of_species)) #num of NC with all negative eigenvalues
        
        # % of performance difference (NC Target Species>0 when param omitted/NC Target Species>0 when no params omitted)*100
        #NOT WORKING
        
        ParameterSensitivity[p+1,6]<-as.numeric(ParameterSensitivity[p+1,4])/as.numeric(ParameterSensitivity[1,4])*100 #percent performance difference with parameter omitted
        
        # num of NC [0 in single species inoc] -> [>0 in probiotic support community inoc]
        #1. [Single Killer Fin States - species 6] AND [final States Killer Comm - species 6] into same matrix, two columns
        perform_metrics_table<-matrix(nrow=num_of_comms, ncol=5)
        perform_metrics_table[,1]<-final_states_KillerComm[,6]
        perform_metrics_table[,2]<-SingleKillerFinStates[,6]
        #2. column 3, give a 1 if single killer Fin state = 0
        perform_metrics_table[,3]<-ifelse(perform_metrics_table[,2]==0,1,0)
        #3. column 4, give a 1 if [killercomm-singleinoc]>0
        perform_metrics_table[,4]<-ifelse(perform_metrics_table[,1]>perform_metrics_table[,2],1,0)
        #4. column 5, add two columns
        perform_metrics_table[,5]<-perform_metrics_table[,3]+perform_metrics_table[,4]
        #5. NC of [0] -> [>0] is sums(column 5 = 2) - number of rows whose column #5 is = 2
        NCs_zero_to_pos<-sum(perform_metrics_table[,5]==2)
        #6. Put this # of NC inot ParameterSensitivity[p+1,7]
        ParameterSensitivity[p+1,8]<-NCs_zero_to_pos
        pct_diff_zeroToPos<-(NCs_zero_to_pos/as.numeric(ParameterSensitivity[1,8]))*100 #calculating % difference in [0 in single species inoc] -> [>0 in probiotic support community inoc] with parameter omitted
        ParameterSensitivity[p+1,9]<-pct_diff_zeroToPos
        
        # of NC [>0 in single species inoc] -> [>>0 in probiotic support community inoc] 
        #1. [Single Killer Fin States - species 6] AND [final States Killer Comm - species 6] into same matrix, two columns
        perform_metrics_table2<-matrix(nrow=num_of_comms, ncol=5)
        perform_metrics_table2[,1]<-final_states_KillerComm[,6] #final state of target species from PSC
        perform_metrics_table2[,2]<-SingleKillerFinStates[,6] #final states of target species from single target species inoculation
        #2. column 3, give a 1 if single killer Fin state > 0
        perform_metrics_table2[,3]<-ifelse(perform_metrics_table2[,2]>0,1,0)
        #3. column 4, give a 1 if [killercomm-singleinoc]>0
        perform_metrics_table2[,4]<-ifelse(perform_metrics_table2[,1]>perform_metrics_table2[,2],1,0)
        #4. column 5, add two columns
        perform_metrics_table2[,5]<-perform_metrics_table2[,3]+perform_metrics_table2[,4]
        #5. NC of [0] -> [>0] is sums(column 5 = 2) - number of rows whose column #5 is = 2
        NCs_pos_to_greater<-sum(perform_metrics_table2[,5]==2)
        #6. Put this # of NC inot ParameterSensitivity[p+1,8]
        ParameterSensitivity[p+1,10]<-NCs_pos_to_greater
        pct_diff_posToGreater<-(NCs_pos_to_greater/as.numeric(ParameterSensitivity[1,10]))*100 #calculating % difference in [>0 in single species inoc] -> [>>0 in probiotic support community inoc] with parameter omitted
        ParameterSensitivity[p+1,11]<-pct_diff_posToGreater
        
        rm(perform_metrics_table)
        rm(perform_metrics_table2)
        
        p = p+1
        print(p)
        n = 1
      }
      
      ParameterSensitivity[2,7]<-paste("a17")
      ParameterSensitivity[3,7]<-paste("a18")
      ParameterSensitivity[4,7]<-paste("a19")
      ParameterSensitivity[5,7]<-paste("a110")
      ParameterSensitivity[6,7]<-paste("a1_11")
      
      ParameterSensitivity[7,7]<-paste("a27")
      ParameterSensitivity[8,7]<-paste("a28")
      ParameterSensitivity[9,7]<-paste("a29")
      ParameterSensitivity[10,7]<-paste("a210")
      ParameterSensitivity[11,7]<-paste("a211")
      
      ParameterSensitivity[12,7]<-paste("a37")
      ParameterSensitivity[13,7]<-paste("a38")
      ParameterSensitivity[14,7]<-paste("a39")
      ParameterSensitivity[15,7]<-paste("a310")
      ParameterSensitivity[16,7]<-paste("a311")
      
      ParameterSensitivity[17,7]<-paste("a47")
      ParameterSensitivity[18,7]<-paste("a48")
      ParameterSensitivity[19,7]<-paste("a49")
      ParameterSensitivity[20,7]<-paste("a410")
      ParameterSensitivity[21,7]<-paste("a411")
      
      ParameterSensitivity[22,7]<-paste("a57")
      ParameterSensitivity[23,7]<-paste("a58")
      ParameterSensitivity[24,7]<-paste("a59")
      ParameterSensitivity[25,7]<-paste("a510")
      ParameterSensitivity[26,7]<-paste("a511")
      
      ParameterSensitivity[27,7]<-paste("a67")
      ParameterSensitivity[28,7]<-paste("a68")
      ParameterSensitivity[29,7]<-paste("a69")
      ParameterSensitivity[30,7]<-paste("a610")
      ParameterSensitivity[31,7]<-paste("a611")
      
      ParameterSensitivity[32,7]<-paste("a71")
      ParameterSensitivity[33,7]<-paste("a72")
      ParameterSensitivity[34,7]<-paste("a73")
      ParameterSensitivity[35,7]<-paste("a74")
      ParameterSensitivity[36,7]<-paste("a75")
      ParameterSensitivity[37,7]<-paste("a76")
      ParameterSensitivity[38,7]<-paste("a78")
      ParameterSensitivity[39,7]<-paste("a79")
      ParameterSensitivity[40,7]<-paste("a710")
      ParameterSensitivity[41,7]<-paste("a711")
      
      ParameterSensitivity[42,7]<-paste("a81")
      ParameterSensitivity[43,7]<-paste("a82")
      ParameterSensitivity[44,7]<-paste("a83")
      ParameterSensitivity[45,7]<-paste("a84")
      ParameterSensitivity[46,7]<-paste("a85")
      ParameterSensitivity[47,7]<-paste("a86")
      ParameterSensitivity[48,7]<-paste("a87")
      ParameterSensitivity[49,7]<-paste("a89")
      ParameterSensitivity[50,7]<-paste("a810")
      ParameterSensitivity[51,7]<-paste("a811")
      
      ParameterSensitivity[52,7]<-paste("a91")
      ParameterSensitivity[53,7]<-paste("a92")
      ParameterSensitivity[54,7]<-paste("a93")
      ParameterSensitivity[55,7]<-paste("a94")
      ParameterSensitivity[56,7]<-paste("a95")
      ParameterSensitivity[57,7]<-paste("a96")
      ParameterSensitivity[58,7]<-paste("a97")
      ParameterSensitivity[59,7]<-paste("a98")
      ParameterSensitivity[60,7]<-paste("a910")
      ParameterSensitivity[61,7]<-paste("a911")
      
      ParameterSensitivity[62,7]<-paste("a101")
      ParameterSensitivity[63,7]<-paste("a102")
      ParameterSensitivity[64,7]<-paste("a103")
      ParameterSensitivity[65,7]<-paste("a104")
      ParameterSensitivity[66,7]<-paste("a105")
      ParameterSensitivity[67,7]<-paste("a106")
      ParameterSensitivity[68,7]<-paste("a107")
      ParameterSensitivity[69,7]<-paste("a108")
      ParameterSensitivity[70,7]<-paste("a109")
      ParameterSensitivity[71,7]<-paste("a1011")
      
      ParameterSensitivity[72,7]<-paste("a11_1")
      ParameterSensitivity[73,7]<-paste("a112")
      ParameterSensitivity[74,7]<-paste("a113")
      ParameterSensitivity[75,7]<-paste("a114")
      ParameterSensitivity[76,7]<-paste("a115")
      ParameterSensitivity[77,7]<-paste("a116")
      ParameterSensitivity[78,7]<-paste("a117")
      ParameterSensitivity[79,7]<-paste("a118")
      ParameterSensitivity[80,7]<-paste("a119")
      ParameterSensitivity[81,7]<-paste("a1110")
      
      #This command likely needs to be updated......see what the table looks like
      #update file write location - into set file directory
      #DO write column names
      writeFilePath<-paste("~/Documents/SynCom_Modeling/",file_number,"/ParameterSensitivity_PSC_",PSC_number,".csv",sep="")
      write.table(ParameterSensitivity, file = writeFilePath, sep = ",",row.names=FALSE)
