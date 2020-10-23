rm(lists=ls())
file_number<-4
dataframeHeader<- read.csv("~/Documents/Arkin_Lab/CUBES/Experiments_Data_Planning_Priortizing/Microbiome Engineering/community modeling/R package_our code documentation/13apr20_PSC_parameter_files_headers.csv", header=FALSE) #read in file containing header names
dataframeHeadTrans<-t(dataframeHeader) #transpose header file
dataframeHeadTrans<-dataframeHeadTrans[,2:133]
PSCParams<-matrix(ncol=132,nrow=2000)
colnames(PSCParams)<-(dataframeHeadTrans[2,]) #generating first row of dataframeImport matrix
PSCParams<-as.data.frame(PSCParams)
for(i in 1:2000){
  #make filepath with icsvFileName <- paste("KillerComm_",p,"_stode_final_parameters.csv",sep="")
  csvFileName <- paste("KillerComm_",i,"_stode_final_parameters.csv",sep="")
  csvFilePath <-paste("~/Documents/SynCom_Modeling/",file_number,"/KillerComm_stode_community_parameters/",csvFileName,sep="") 
  #importParams<-read.csv(file = csvFilePath, sep = ",", header=F, col.names = F, row.names = F) #read in PSC full params list into importParams
  importParams <- read.csv(csvFilePath, header=FALSE, row.names=NULL, stringsAsFactors=FALSE)
  PSCParams[i,]<-importParams[1,] #copy first row of imported params file to dataFrameImport
  #rm(importParams)
  print(i)
}

PSCParamsMat<-as.matrix(PSCParams, colnames=FALSE)
ParamHeadings<-read.csv("~/Desktop/reconfigParamHeadings.csv")
PSCParams_Reconfig<-matrix(nrow=2000, ncol=80)
colnames(PSCParams_Reconfig)<-colnames(ParamHeadings)

j<-1
k<-1
for (j in 1:80){
  for (k in 1:132){
    if (colnames(PSCParams_Reconfig)[j]==colnames(PSCParamsMat)[k]){
      PSCParams_Reconfig[,j]<-PSCParamsMat[,k]
      #put column k from PSC params as column j in PSCParams_reconfig
  
    }
  
  }
}

write.csv(PSCParams_Reconfig,"~/Desktop/PSCParams_eigenstable_reconfigured.csv") #write table to desktop for import into libreoffice calc
    