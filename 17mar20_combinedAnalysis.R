library(rootSolve)
 
# This script is Used to analyze native communities, Single Killer Communities, and KillerComms

rm(list=ls())
#Fill These Out
num_of_KillerComm<-2000
num_of_comms<-1000
file_number<-3 #file number name where KillerComm and Single Killer files are stored
num_of_species<-11 #total number of species in final combined KillerComm

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

#header names - add to table
#left off here - need to add two columns of 0->something and >0 -> >>0 to ave, min and max
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

#add headers to combinedMetrics.csv
