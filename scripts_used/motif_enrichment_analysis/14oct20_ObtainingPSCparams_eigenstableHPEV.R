rm(lists=ls())

PSCParamsMat<-read.csv("~/Documents/Arkin_Lab/CUBES/Experiments_Data_Planning_Priortizing/Microbiome Engineering/community modeling/functional evolution_enrichment analysis/eigenstable reconfiguring parameter matricies/14oct20_HP_EV_eigentstable evolution_params NOT reconfigured.csv")
PSCParamsMat<-as.matrix(PSCParamsMat)
ParamHeadings<-read.csv("~/Documents/Arkin_Lab/CUBES/Experiments_Data_Planning_Priortizing/Microbiome Engineering/community modeling/functional evolution_enrichment analysis/eigenstable reconfiguring parameter matricies/reconfigParamHeadings.csv")
PSCParams_Reconfig<-matrix(nrow=12, ncol=80)
colnames(PSCParams_Reconfig)<-colnames(ParamHeadings)
rownames(PSCParams_Reconfig<-rownames(PSCParamsMat))
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
 
  write.csv(PSCParams_Reconfig,"~/Desktop/PSCParams_eigenstable_HP_EV_reconfigured.csv") #write table to desktop for import into libreoffice calc
  