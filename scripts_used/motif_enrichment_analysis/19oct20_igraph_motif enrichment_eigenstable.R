install.packages("igraph")
library(igraph)
#sometimes after installing and loading packages/libraries for igraph, the igraph functions give errors
#reinstalling package and reloading library a second time seems to fix this

rm(list=ls())
file_number<-4 #4 for eigenstable initial PSC's
num_of_KillerComm <- 2000
num_of_comms <-1000
num_of_species<-11
final_parameters_KillerComm <- matrix(ncol = (num_of_species*num_of_species+num_of_species), nrow = num_of_comms)
motif_counts_3node<-matrix(ncol=16,nrow=num_of_KillerComm)
motif_counts_4node<-matrix(ncol=218,nrow=num_of_KillerComm)
motif_counts_3node_5pct<-matrix(ncol=16,nrow=num_of_KillerComm)
motif_counts_4node_5pct<-matrix(ncol=218,nrow=num_of_KillerComm)
motif_counts_3node_10pct<-matrix(ncol=16,nrow=num_of_KillerComm)
motif_counts_4node_10pct<-matrix(ncol=218,nrow=num_of_KillerComm)
motif_counts_3node_15pct<-matrix(ncol=16,nrow=num_of_KillerComm)
motif_counts_4node_15pct<-matrix(ncol=218,nrow=num_of_KillerComm)
motif_countsHP_3node<-matrix(ncol=16,nrow=5)
motif_countsHP_4node<-matrix(ncol=218,nrow=5)
motif_countsHP_3node_5pct<-matrix(ncol=16,nrow=5)
motif_countsHP_4node_5pct<-matrix(ncol=218,nrow=5)
motif_countsHP_3node_10pct<-matrix(ncol=16,nrow=5)
motif_countsHP_4node_10pct<-matrix(ncol=218,nrow=5)
motif_countsHP_3node_15pct<-matrix(ncol=16,nrow=5)
motif_countsHP_4node_15pct<-matrix(ncol=218,nrow=5)
motif_countsEV_3node<-matrix(ncol=16,nrow=5)
motif_countsEV_4node<-matrix(ncol=218,nrow=5)
motif_countsEV_3node_5pct<-matrix(ncol=16,nrow=5)
motif_countsEV_4node_5pct<-matrix(ncol=218,nrow=5)
motif_countsEV_3node_10pct<-matrix(ncol=16,nrow=5)
motif_countsEV_4node_10pct<-matrix(ncol=218,nrow=5)
motif_countsEV_3node_15pct<-matrix(ncol=16,nrow=5)
motif_countsEV_4node_15pct<-matrix(ncol=218,nrow=5)


#Whole PSC population from stode native comms - motif counts, isomorphism class means and stdevs
for (p in 1:num_of_KillerComm){     
csvFileName <- paste("KillerComm_",p,"_stode_final_parameters.csv",sep="")
csvFilePath <-paste("~/Documents/SynCom_Modeling/",file_number,"/KillerComm_stode_community_parameters/",csvFileName,sep="") 
final_parameters_KillerComm<-as.matrix(read.csv(csvFilePath, sep = ",", header = FALSE))
#make adjacency matrix
adj_mat_T<-matrix(ncol=11,nrow=11) #initializes transpose of adjacency matrix - this matrix will be filled in from parameter files
adj_mat_T[1:6,7:11]<-rbind(final_parameters_KillerComm[1,18:22], final_parameters_KillerComm[1,29:33], final_parameters_KillerComm[1,40:44], final_parameters_KillerComm[1,51:55], final_parameters_KillerComm[1,62:66], final_parameters_KillerComm[1,73:77])#filling in rows 1-6 of adj_mat_T
adj_mat_T[7:11,1:11]<-rbind(final_parameters_KillerComm[1,78:88], final_parameters_KillerComm[1,89:99], final_parameters_KillerComm[1,100:110], final_parameters_KillerComm[1,111:121], final_parameters_KillerComm[1,122:132])#filling in rows 7-11 of adj_mat_T
adj_mat_T[1:6,1:6]<-c(0,0,0,0,0,0)#filling in zeros on adj_matrix_T
adj_mat<-t(adj_mat_T) #transpose adjacency matrix - have to do this for igraph
ZeroPctGraph<-graph_from_adjacency_matrix(adj_mat, mode=c("directed"), weighted=TRUE) #make graph from adjacency matrix
motif_counts_3node[p,1:16]<-motifs(ZeroPctGraph,size=3) #count 3 node motifs, and output count vector to results table
motif_counts_4node[p,1:218]<-motifs(ZeroPctGraph,size=4) #count 4 node motifs, and output count vector to results table
#NA values in vector refer to isomorphism classes that only include 2 nodes, not 3 nodes, thus the NA values

#below abs(5%) of max(adj_mat) entry, set to zero

if (abs(min(adj_mat))>max(adj_mat)){
v<-0.05*(abs(min(adj_mat)))
} else{
v<-0.05*(max(adj_mat))
}

adj_mat_5<-matrix(ncol=11,nrow=11) #initializes transpose of adjacency matrix for 5% set to zero - this matrix will be filled in by 5% calculations
i<-1
j<-1
#set all values < abs(5%) of max(adj_mat) to zero
for (i in 1:11){
  for (j in 1:11){
    if((adj_mat[i,j]<v) & (adj_mat[i,j]>(v*-1))){
      adj_mat_5[i,j]<-0}
    else{
      adj_mat_5[i,j]<-adj_mat[i,j]  
    }
  }
}
FivepctGraph<-graph_from_adjacency_matrix(adj_mat_5, mode=c("directed"), weighted=TRUE)
motif_counts_3node_5pct[p,1:16]<-motifs(FivepctGraph,size=3) #count 3 node motifs, and output count vector to results table
motif_counts_4node_5pct[p,1:218]<-motifs(FivepctGraph,size=4) #count 4 node motifs, and output count vector to results table

#below abs(10%) set to zero
if (abs(min(adj_mat))>max(adj_mat)){ #if statement to determine if max(adj_mat)<min(abs(adj_mat))
  v<-0.1*(abs(min(adj_mat)))
} else{
  v<-0.1*(max(adj_mat))
}
adj_mat_10<-matrix(ncol=11,nrow=11) #initializes transpose of adjacency matrix for 10% set to zero - this matrix will be filled in by 10% comparison
i<-1
j<-1
#set all values < abs(10%) of max(adj_mat) to zero
for (i in 1:11){
  for (j in 1:11){
    if((adj_mat[i,j]<v) & (adj_mat[i,j]>(v*-1))){
      adj_mat_10[i,j]<-0}
    else{
      adj_mat_10[i,j]<-adj_mat[i,j]  
    }
  }
}
TenPctGraph<-graph_from_adjacency_matrix(adj_mat_10, mode=c("directed"), weighted=TRUE)
motif_counts_3node_10pct[p,1:16]<-motifs(TenPctGraph,size=3) #count 3 node motifs, and output count vector to results table
motif_counts_4node_10pct[p,1:218]<-motifs(TenPctGraph,size=4) #count 4 node motifs, and output count vector to results table

#below abs(15%) set to zero
if (abs(min(adj_mat))>max(adj_mat)){
  v<-0.15*(abs(min(adj_mat)))
} else{
  v<-0.15*(max(adj_mat))
}
adj_mat_15<-matrix(ncol=11,nrow=11) #initializes transpose of adjacency matrix for 15% set to zero - this matrix will be filled in by 10% comparison
i<-1
j<-1
#set all values < abs(10%) of max(adj_mat) to zero
for (i in 1:11){
  for (j in 1:11){
    if((adj_mat[i,j]<v) & (adj_mat[i,j]>(v*-1))){
      adj_mat_15[i,j]<-0}
    else{
      adj_mat_15[i,j]<-adj_mat[i,j]  
    }
  }
}
FifteenPctGraph<-graph_from_adjacency_matrix(adj_mat_15, mode=c("directed"), weighted=TRUE)
motif_counts_3node_15pct[p,1:16]<-motifs(FifteenPctGraph,size=3) #count 3 node motifs, and output count vector to results table
motif_counts_4node_15pct[p,1:218]<-motifs(FifteenPctGraph,size=4) #count 4 node motifs, and output count vector to results table
print(p)
} #end of p for loop for all 2000 PSC's

#compiles 3-node averages for 0%, 5%, 10%, 15% into matrix (a matrix, rows are isomorphism class, columns are 0%, 5%, 10%, 15% thresholds)
ThreeNodeMeansPop<-cbind(colMeans(motif_counts_3node, na.rm=TRUE), colMeans(motif_counts_3node_5pct, na.rm=TRUE), colMeans(motif_counts_3node_10pct, na.rm=TRUE), colMeans(motif_counts_3node_15pct, na.rm=TRUE)) #gives column means for each isomorphism class
threenodemeancolnamesPop<-c('Pop Mean 0%','Pop Mean 5%','Pop Mean 10%','Pop Mean 15%')
colnames(ThreeNodeMeansPop)<-threenodemeancolnamesPop #adds colnames to ThreeNodeMeans

#compiles 3-node standard deviations for 0%, 5%, 10%, 15% into matrix (a matrix, rows are isomorphism class, columns are 0%, 5%, 10%, 15% thresholds)
ThreeNodeSDsPop<-cbind(apply(motif_counts_3node,2,sd, na.rm=TRUE),apply(motif_counts_3node_5pct,2,sd, na.rm=TRUE), apply(motif_counts_3node_10pct,2,sd, na.rm=TRUE), apply(motif_counts_3node_15pct,2,sd, na.rm=TRUE)) #gives column SD's for each isomorphism class
threenodeSDcolnamesPop<-c('Pop SD 0%','Pop SD 5%','Pop SD 10%','Pop SD 15%')
colnames(ThreeNodeSDsPop)<-threenodeSDcolnamesPop #add colnames to ThreeNodeSDs 

#compiles 4-node averages for 0%, 5%, 10%, 15% into matrix (a matrix, rows are isomorphism class, columns are 0%, 5%, 10%, 15% thresholds)
FourNodeMeansPop<-cbind(colMeans(motif_counts_4node, na.rm=TRUE), colMeans(motif_counts_4node_5pct, na.rm=TRUE), colMeans(motif_counts_4node_10pct, na.rm=TRUE), colMeans(motif_counts_4node_15pct, na.rm=TRUE)) #gives column means for each isomorphism class
fournodemeancolnamesPop<-c('Pop Mean 0%','Pop Mean 5%','Pop Mean 10%','Pop Mean 15%')
colnames(FourNodeMeansPop)<-fournodemeancolnamesPop #adds colnames to FourNodeMeans

#compiles 4-node standard deviations for 0%, 5%, 10%, 15% into matrix (a matrix, rows are isomorphism class, columns are 0%, 5%, 10%, 15% thresholds)
FourNodeSDsPop<-cbind(apply(motif_counts_4node,2,sd, na.rm=TRUE),apply(motif_counts_4node_5pct,2,sd, na.rm=TRUE), apply(motif_counts_4node_10pct,2,sd, na.rm=TRUE), apply(motif_counts_4node_15pct,2,sd, na.rm=TRUE)) #gives column SD's for each isomorphism class
fournodeSDcolnamesPop<-c('Pop SD 0%','Pop SD 5%','Pop SD 10%','Pop SD 15%')
colnames(FourNodeSDsPop)<-fournodeSDcolnamesPop #add colnames to FourNodeSDsA
  
#Top 5 high-performing PSCs (HP) from stode native comms - 3 and 4 node motif counts, means, stdevs
p<-1
for (p in 1:5){ #reading in top 5 high performing communities from stode PSCs
  csvFileName <- paste("KillerComm_",p,"_stode_final_parameters.csv",sep="")
  csvFilePath <-paste("~/Documents/SynCom_Modeling/param_files_enrichment/stode/high_perf/",csvFileName,sep="") 
  final_parameters_KillerComm<-as.matrix(read.csv(csvFilePath, sep = ",", header = FALSE))
  adj_mat_T[1:6,7:11]<-rbind(final_parameters_KillerComm[1,18:22], final_parameters_KillerComm[1,29:33], final_parameters_KillerComm[1,40:44], final_parameters_KillerComm[1,51:55], final_parameters_KillerComm[1,62:66], final_parameters_KillerComm[1,73:77])#filling in rows 1-6 of adj_mat_T
  adj_mat_T[7:11,1:11]<-rbind(final_parameters_KillerComm[1,78:88], final_parameters_KillerComm[1,89:99], final_parameters_KillerComm[1,100:110], final_parameters_KillerComm[1,111:121], final_parameters_KillerComm[1,122:132])#filling in rows 7-11 of adj_mat_T
  adj_mat_T[1:6,1:6]<-c(0,0,0,0,0,0)#filling in zeros on adj_matrix_T
  adj_mat<-t(adj_mat_T) #transpose adjacency matrix - have to do this - igraph indicies convention
  ZeroPctGraphHP<-graph_from_adjacency_matrix(adj_mat, mode=c("directed"), weighted=TRUE) #make graph from adjacency matrix
  motif_countsHP_3node[p,1:16]<-motifs(ZeroPctGraphHP,size=3) #count 3 node motifs, and output count vector to results table
  motif_countsHP_4node[p,1:218]<-motifs(ZeroPctGraphHP,size=4) #count 4 node motifs, and output count vector to results table

  #below abs(5%) of max(adj_mat) entry, set to zero
  
  if (abs(min(adj_mat))>max(adj_mat)){
    v<-0.05*(abs(min(adj_mat)))
  } else{
    v<-0.05*(max(adj_mat))
  }
  i<-1
  j<-1
  #set all values < abs(5%) of max(adj_mat) to zero
  for (i in 1:11){
    for (j in 1:11){
      if((adj_mat[i,j]<v) & (adj_mat[i,j]>(v*-1))){
        adj_mat_5[i,j]<-0}
      else{
        adj_mat_5[i,j]<-adj_mat[i,j]  
      }
    }
  }
  FivepctGraphHP<-graph_from_adjacency_matrix(adj_mat_5, mode=c("directed"), weighted=TRUE)
  motif_countsHP_3node_5pct[p,1:16]<-motifs(FivepctGraphHP,size=3) #count 3 node motifs, and output count vector to results table
  motif_countsHP_4node_5pct[p,1:218]<-motifs(FivepctGraphHP,size=4) #count 4 node motifs, and output count vector to results table
  
  #below abs(10%) set to zero
  if (abs(min(adj_mat))>max(adj_mat)){
    v<-0.1*(abs(min(adj_mat)))
  } else{
    v<-0.1*(max(adj_mat))
  }
  i<-1
  j<-1
  #set all values < abs(10%) of max(adj_mat) to zero
  for (i in 1:11){
    for (j in 1:11){
      if((adj_mat[i,j]<v) & (adj_mat[i,j]>(v*-1))){
        adj_mat_10[i,j]<-0}
      else{
        adj_mat_10[i,j]<-adj_mat[i,j]  
      }
    }
  }
  TenPctGraphHP<-graph_from_adjacency_matrix(adj_mat_10, mode=c("directed"), weighted=TRUE)
  motif_countsHP_3node_10pct[p,1:16]<-motifs(TenPctGraphHP,size=3) #count 3 node motifs, and output count vector to results table
  motif_countsHP_4node_10pct[p,1:218]<-motifs(TenPctGraphHP,size=4) #count 4 node motifs, and output count vector to results table
  
  #below abs(15%) set to zero
  if (abs(min(adj_mat))>max(adj_mat)){
    v<-0.15*(abs(min(adj_mat)))
  } else{
    v<-0.15*(max(adj_mat))
  }
  i<-1
  j<-1
  #set all values < abs(10%) of max(adj_mat) to zero
  for (i in 1:11){
    for (j in 1:11){
      if((adj_mat[i,j]<v) & (adj_mat[i,j]>(v*-1))){
        adj_mat_15[i,j]<-0}
      else{
        adj_mat_15[i,j]<-adj_mat[i,j]  
      }
    }
  }
  FifteenPctGraphHP<-graph_from_adjacency_matrix(adj_mat_15, mode=c("directed"), weighted=TRUE)
  motif_countsHP_3node_15pct[p,1:16]<-motifs(FifteenPctGraphHP,size=3) #count 3 node motifs, and output count vector to results table
  motif_countsHP_4node_15pct[p,1:218]<-motifs(FifteenPctGraphHP,size=4) #count 4 node motifs, and output count vector to results table
  print(p)
} #end of p for loop for top 5 high-performing PSC's

#compiles 3-node averages for 0%, 5%, 10%, 15% into matrix (a matrix, rows are isomorphism class, columns are 0%, 5%, 10%, 15% thresholds)
ThreeNodeMeansHP<-cbind(colMeans(motif_countsHP_3node, na.rm=TRUE), colMeans(motif_countsHP_3node_5pct, na.rm=TRUE), colMeans(motif_countsHP_3node_10pct, na.rm=TRUE), colMeans(motif_countsHP_3node_15pct, na.rm=TRUE)) #gives column means for each isomorphism class
threenodemeancolnamesHP<-c('HP Mean 0%','HP Mean 5%','HP Mean 10%','HP Mean 15%')
colnames(ThreeNodeMeansHP)<-threenodemeancolnamesHP #adds colnames to ThreeNodeMeansHP

#compiles 3-node standard deviations for 0%, 5%, 10%, 15% into matrix (a matrix, rows are isomorphism class, columns are 0%, 5%, 10%, 15% thresholds)
ThreeNodeSDsHP<-cbind(apply(motif_countsHP_3node,2,sd, na.rm=TRUE),apply(motif_countsHP_3node_5pct,2,sd, na.rm=TRUE), apply(motif_countsHP_3node_10pct,2,sd, na.rm=TRUE), apply(motif_countsHP_3node_15pct,2,sd, na.rm=TRUE)) #gives column SD's for each isomorphism class
threenodeSDcolnamesHP<-c('HP SD 0%','HP SD 5%','HP SD 10%','HP SD 15%')
colnames(ThreeNodeSDsHP)<-threenodeSDcolnamesHP #add colnames to ThreeNodeSDsHP 

#compiles 4-node averages for 0%, 5%, 10%, 15% into matrix (a matrix, rows are isomorphism class, columns are 0%, 5%, 10%, 15% thresholds)
FourNodeMeansHP<-cbind(colMeans(motif_countsHP_4node, na.rm=TRUE), colMeans(motif_countsHP_4node_5pct, na.rm=TRUE), colMeans(motif_countsHP_4node_10pct, na.rm=TRUE), colMeans(motif_countsHP_4node_15pct, na.rm=TRUE)) #gives column means for each isomorphism class
fournodemeancolnamesHP<-c('HP Mean 0%','HP Mean 5%','HP Mean 10%','HP Mean 15%')
colnames(FourNodeMeansHP)<-fournodemeancolnamesHP #adds colnames to FourNodeMeansHP

#compiles 4-node standard deviations for 0%, 5%, 10%, 15% into matrix (a matrix, rows are isomorphism class, columns are 0%, 5%, 10%, 15% thresholds)
FourNodeSDsHP<-cbind(apply(motif_countsHP_4node,2,sd, na.rm=TRUE),apply(motif_countsHP_4node_5pct,2,sd, na.rm=TRUE), apply(motif_countsHP_4node_10pct,2,sd, na.rm=TRUE), apply(motif_countsHP_4node_15pct,2,sd, na.rm=TRUE)) #gives column SD's for each isomorphism class
fournodeSDcolnamesHP<-c('HP SD 0%','HP SD 5%','HP SD 10%','HP SD 15%')
colnames(FourNodeSDsHP)<-fournodeSDcolnamesHP #add colnames to FourNodeSDsHP

#Top 5 evolved PSCs (HP) from stode native comms - 3 and 4 node motif counts, means, stdevs
p<-1
for (p in 1:5){ #reading in top 5 high performing communities from stode PSCs
  csvFileName <- paste("KillerComm_",p,"_stode_final_parameters.csv",sep="")
  csvFilePath <-paste("~/Documents/SynCom_Modeling/param_files_enrichment/stode/evolved/",csvFileName,sep="") 
  final_parameters_KillerComm<-as.matrix(read.csv(csvFilePath, sep = ",", header = FALSE))
  adj_mat_T[1:6,7:11]<-rbind(final_parameters_KillerComm[1,18:22], final_parameters_KillerComm[1,29:33], final_parameters_KillerComm[1,40:44], final_parameters_KillerComm[1,51:55], final_parameters_KillerComm[1,62:66], final_parameters_KillerComm[1,73:77])#filling in rows 1-6 of adj_mat_T
  adj_mat_T[7:11,1:11]<-rbind(final_parameters_KillerComm[1,78:88], final_parameters_KillerComm[1,89:99], final_parameters_KillerComm[1,100:110], final_parameters_KillerComm[1,111:121], final_parameters_KillerComm[1,122:132])#filling in rows 7-11 of adj_mat_T
  adj_mat_T[1:6,1:6]<-c(0,0,0,0,0,0)#filling in zeros on adj_matrix_T
  adj_mat<-t(adj_mat_T) #transpose adjacency matrix - have to do this - igraph indicies convention
  ZeroPctGraphEV<-graph_from_adjacency_matrix(adj_mat, mode=c("directed"), weighted=TRUE) #make graph from adjacency matrix
  motif_countsEV_3node[p,1:16]<-motifs(ZeroPctGraphEV,size=3) #count 3 node motifs, and output count vector to results table
  motif_countsEV_4node[p,1:218]<-motifs(ZeroPctGraphEV,size=4) #count 4 node motifs, and output count vector to results table
  
  #below abs(5%) of max(adj_mat) entry, set to zero
  if (abs(min(adj_mat))>max(adj_mat)){
    v<-0.05*(abs(min(adj_mat)))
  } else{
    v<-0.05*(max(adj_mat))
  }
  i<-1
  j<-1
  #set all values < abs(5%) of max(adj_mat) to zero
  for (i in 1:11){
    for (j in 1:11){
      if((adj_mat[i,j]<v) & (adj_mat[i,j]>(v*-1))){
        adj_mat_5[i,j]<-0}
      else{
        adj_mat_5[i,j]<-adj_mat[i,j]  
      }
    }
  }
  FivepctGraphEV<-graph_from_adjacency_matrix(adj_mat_5, mode=c("directed"), weighted=TRUE)
  motif_countsEV_3node_5pct[p,1:16]<-motifs(FivepctGraphEV,size=3) #count 3 node motifs, and output count vector to results table
  motif_countsEV_4node_5pct[p,1:218]<-motifs(FivepctGraphEV,size=4) #count 4 node motifs, and output count vector to results table
  
  #below abs(10%) set to zero
  if (abs(min(adj_mat))>max(adj_mat)){
    v<-0.1*(abs(min(adj_mat)))
  } else{
    v<-0.1*(max(adj_mat))
  }
  i<-1
  j<-1
  #set all values < abs(10%) of max(adj_mat) to zero
  for (i in 1:11){
    for (j in 1:11){
      if((adj_mat[i,j]<v) & (adj_mat[i,j]>(v*-1))){
        adj_mat_10[i,j]<-0}
      else{
        adj_mat_10[i,j]<-adj_mat[i,j]  
      }
    }
  }
  TenPctGraphEV<-graph_from_adjacency_matrix(adj_mat_10, mode=c("directed"), weighted=TRUE)
  motif_countsEV_3node_10pct[p,1:16]<-motifs(TenPctGraphEV,size=3) #count 3 node motifs, and output count vector to results table
  motif_countsEV_4node_10pct[p,1:218]<-motifs(TenPctGraphEV,size=4) #count 4 node motifs, and output count vector to results table
  
  #below abs(15%) set to zero
  if (abs(min(adj_mat))>max(adj_mat)){
    v<-0.15*(abs(min(adj_mat)))
  } else{
    v<-0.15*(max(adj_mat))
  }
  i<-1
  j<-1
  #set all values < abs(10%) of max(adj_mat) to zero
  for (i in 1:11){
    for (j in 1:11){
      if((adj_mat[i,j]<v) & (adj_mat[i,j]>(v*-1))){
        adj_mat_15[i,j]<-0}
      else{
        adj_mat_15[i,j]<-adj_mat[i,j]  
      }
    }
  }
  FifteenPctGraphEV<-graph_from_adjacency_matrix(adj_mat_15, mode=c("directed"), weighted=TRUE)
  motif_countsEV_3node_15pct[p,1:16]<-motifs(FifteenPctGraphEV,size=3) #count 3 node motifs, and output count vector to results table
  motif_countsEV_4node_15pct[p,1:218]<-motifs(FifteenPctGraphEV,size=4) #count 4 node motifs, and output count vector to results table
  print(p)
}#end of p for loop for top 5 evolved PSC's

#compiles 3-node averages for 0%, 5%, 10%, 15% into matrix (a matrix, rows are isomorphism class, columns are 0%, 5%, 10%, 15% thresholds)
ThreeNodeMeansEV<-cbind(colMeans(motif_countsEV_3node, na.rm=TRUE), colMeans(motif_countsEV_3node_5pct, na.rm=TRUE), colMeans(motif_countsEV_3node_10pct, na.rm=TRUE), colMeans(motif_countsEV_3node_15pct, na.rm=TRUE)) #gives column means for each isomorphism class
threenodemeancolnamesEV<-c('EV Mean 0%','EV Mean 5%','EV Mean 10%','EV Mean 15%')
colnames(ThreeNodeMeansEV)<-threenodemeancolnamesEV #adds colnames to ThreeNodeMeansEV

#compiles 3-node standard deviations for 0%, 5%, 10%, 15% into matrix (a matrix, rows are isomorphism class, columns are 0%, 5%, 10%, 15% thresholds)
ThreeNodeSDsEV<-cbind(apply(motif_countsEV_3node,2,sd, na.rm=TRUE),apply(motif_countsEV_3node_5pct,2,sd, na.rm=TRUE), apply(motif_countsEV_3node_10pct,2,sd, na.rm=TRUE), apply(motif_countsEV_3node_15pct,2,sd, na.rm=TRUE)) #gives column SD's for each isomorphism class
threenodeSDcolnamesEV<-c('EV SD 0%','EV SD 5%','EV SD 10%','EV SD 15%')
colnames(ThreeNodeSDsEV)<-threenodeSDcolnamesEV #add colnames to ThreeNodeSDsEV

#compiles 4-node averages for 0%, 5%, 10%, 15% into matrix (a matrix, rows are isomorphism class, columns are 0%, 5%, 10%, 15% thresholds)
FourNodeMeansEV<-cbind(colMeans(motif_countsEV_4node, na.rm=TRUE), colMeans(motif_countsEV_4node_5pct, na.rm=TRUE), colMeans(motif_countsEV_4node_10pct, na.rm=TRUE), colMeans(motif_countsEV_4node_15pct, na.rm=TRUE)) #gives column means for each isomorphism class
fournodemeancolnamesEV<-c('EV Mean 0%','EV Mean 5%','EV Mean 10%','EV Mean 15%')
colnames(FourNodeMeansEV)<-fournodemeancolnamesEV #adds colnames to FourNodeMeansEV

#compiles 4-node standard deviations for 0%, 5%, 10%, 15% into matrix (a matrix, rows are isomorphism class, columns are 0%, 5%, 10%, 15% thresholds)
FourNodeSDsEV<-cbind(apply(motif_countsEV_4node,2,sd, na.rm=TRUE),apply(motif_countsEV_4node_5pct,2,sd, na.rm=TRUE), apply(motif_countsEV_4node_10pct,2,sd, na.rm=TRUE), apply(motif_countsEV_4node_15pct,2,sd, na.rm=TRUE)) #gives column SD's for each isomorphism class
fournodeSDcolnamesEV<-c('EV SD 0%','EV SD 5%','EV SD 10%','EV SD 15%')
colnames(FourNodeSDsEV)<-fournodeSDcolnamesEV #add colnames to FourNodeSDsEV

p<-1

#Calculates z scores and p values for individual enrichment comparisons
ThreeNodeZscoresEV<-(ThreeNodeMeansEV-ThreeNodeMeansPop)/ThreeNodeSDsPop
ThreeNodePvalEV<-pnorm(ThreeNodeZscoresEV)
ThreeNodeZscoresHP<-(ThreeNodeMeansHP-ThreeNodeMeansPop)/ThreeNodeSDsPop
ThreeNodePvalHP<-pnorm(ThreeNodeMeansHP)
FourNodeZscoresEV<-(FourNodeMeansEV-FourNodeMeansPop)/FourNodeSDsPop
FourNodePvalEV<-pnorm(FourNodeZscoresEV)
FourNodeZscoresHP<-(FourNodeMeansHP-FourNodeMeansPop)/FourNodeSDsPop
FourNodePvalHP<-pnorm(FourNodeZscoresHP)

#perform banjamini hochburg adjustment of p values, if many significant - not doing this because none are significant

install.packages("gplots")
library(gplots)
install.packages("RColorBrewer")
library(RColorBrewer)


pdf(file="~/Documents/SynCom_Modeling/param_files_enrichment/stode/FourNodePvalEV_stode.pdf")
heatmap.2(FourNodePvalEV[,2:4],Rowv=NA, Colv = NA, scale = "none", legend, trace="none", density.info = "none", key.xlab="p-value", dendrogram = "none")
dev.off()

pdf(file="~/Documents/SynCom_Modeling/param_files_enrichment/stode/FourNodePvalHP_stode.pdf")
heatmap.2(FourNodePvalHP[,2:4],Rowv=NA, Colv = NA, scale = "none", legend, trace="none", density.info = "none", key.xlab="p-value", dendrogram = "none")
dev.off()

pdf(file="~/Documents/SynCom_Modeling/param_files_enrichment/stode/ThreeNodePvalEV_stode.pdf")
heatmap.2(ThreeNodePvalEV[,2:4],Rowv=NA, Colv = NA, scale = "none", legend, trace="none", density.info = "none", key.xlab="p-value", dendrogram = "none")
dev.off()

pdf(file="~/Documents/SynCom_Modeling/param_files_enrichment/stode/ThreeNodePvalHP_stode.pdf")
heatmap.2(ThreeNodePvalHP[,2:4],Rowv=NA, Colv = NA, scale = "none", legend, trace="none", density.info = "none", key.xlab="p-value", dendrogram = "none")
dev.off()

install.packages("nortest") #installs packages for normality testing, will be used to test normality of count data
library(nortest)
PearsonPval_motif_counts_3node<-matrix(ncol=16,nrow=9)
pvalColNames<-c("3node_5pct","3node_10pct","3node_15pct","3nodeHP_5pct","3nodeHP_10pct","3nodeHP_15pct","3nodeEV_5pct","3nodeEV_10pct","3nodeEV_15pct")
rownames(PearsonPval_motif_counts_3node)<-pvalColNames

#pearson test for normality 3-node motif counts
for (i in 1:16){
  x<-pearson.test(motif_counts_3node_5pct[,i])
  PearsonPval_motif_counts_3node[1,i]<-x$p.value
  y<-pearson.test(motif_counts_3node_10pct[,i])
  PearsonPval_motif_counts_3node[2,i]<-y$p.value
  z<-pearson.test(motif_counts_3node_15pct[,i])
  PearsonPval_motif_counts_3node[3,i]<-z$p.value
  a<-pearson.test(motif_countsHP_3node_5pct[,i])
  PearsonPval_motif_counts_3node[4,i]<-a$p.value
  b<-pearson.test(motif_countsHP_3node_10pct[,i])
  PearsonPval_motif_counts_3node[5,i]<-b$p.value
  c<-pearson.test(motif_countsHP_3node_15pct[,i])
  PearsonPval_motif_counts_3node[6,i]<-c$p.value
  d<-pearson.test(motif_countsEV_3node_5pct[,i])
  PearsonPval_motif_counts_3node[7,i]<-d$p.value
  e<-pearson.test(motif_countsEV_3node_10pct[,i])
  PearsonPval_motif_counts_3node[8,i]<-e$p.value
  f<-pearson.test(motif_countsEV_3node_15pct[,i])
  PearsonPval_motif_counts_3node[9,i]<-f$p.value
}

#pearson test for normality 4-node motif counts
PearsonPval_motif_counts_4node<-matrix(ncol=218,nrow=9)
pvalColNames<-c("4node_5pct","4node_10pct","4node_15pct","4nodeHP_5pct","4nodeHP_10pct","4nodeHP_15pct","4nodeEV_5pct","4nodeEV_10pct","4nodeEV_15pct")
rownames(PearsonPval_motif_counts_4node)<-pvalColNames
i<-1
for (i in 1:218){
  x<-pearson.test(motif_counts_4node_5pct[,i])
  PearsonPval_motif_counts_4node[1,i]<-x$p.value
  y<-pearson.test(motif_counts_4node_5pct[,i])
  PearsonPval_motif_counts_4node[2,i]<-y$p.value
  z<-pearson.test(motif_counts_4node_15pct[,i])
  PearsonPval_motif_counts_4node[3,i]<-z$p.value
  a<-pearson.test(motif_countsHP_4node_5pct[,i])
  PearsonPval_motif_counts_4node[4,i]<-a$p.value
  b<-pearson.test(motif_countsHP_4node_10pct[,i])
  PearsonPval_motif_counts_4node[5,i]<-b$p.value
  c<-pearson.test(motif_countsHP_4node_15pct[,i])
  PearsonPval_motif_counts_4node[6,i]<-c$p.value
  d<-pearson.test(motif_countsEV_4node_5pct[,i])
  PearsonPval_motif_counts_4node[7,i]<-d$p.value
  e<-pearson.test(motif_countsEV_4node_10pct[,i])
  PearsonPval_motif_counts_4node[8,i]<-e$p.value
  f<-pearson.test(motif_countsEV_4node_15pct[,i])
  PearsonPval_motif_counts_4node[9,i]<-f$p.value
}


#qqnorm plots 3node evolved PSCs counts
pdf(file="~/Documents/SynCom_Modeling/param_files_enrichment/qqplots/stode/motif_countsEV_3node_15pct_stode.pdf")
par(mfrow=c(4,4))
for (i in c(3,5:16)){
qqnorm(motif_countsEV_3node_15pct[,i], ylab=paste("IsoClass_",i,"_quantiles",sep=""))
qqline(motif_countsEV_3node_15pct[,i])
}
dev.off()
pdf(file="~/Documents/SynCom_Modeling/param_files_enrichment/qqplots/stode/motif_countsEV_3node_10pct_stode.pdf")
par(mfrow=c(4,4))
for (i in c(3,5:16)){
  qqnorm(motif_countsEV_3node_10pct[,i], ylab=paste("IsoClass_",i,"_quantiles",sep=""))
  qqline(motif_countsEV_3node_10pct[,i])
}
dev.off()
pdf(file="~/Documents/SynCom_Modeling/param_files_enrichment/qqplots/stode/motif_countsEV_3node_5pct_stode.pdf")
par(mfrow=c(4,4))
for (i in c(3,5:16)){
  qqnorm(motif_countsEV_3node_5pct[,i], ylab=paste("IsoClass_",i,"_quantiles",sep=""))
  qqline(motif_countsEV_3node_5pct[,i])
}
dev.off()
#qqnorm plots 3node population counts
pdf(file="~/Documents/SynCom_Modeling/param_files_enrichment/qqplots/stode/motif_countsPopulation_3node_5pct_stode.pdf")
par(mfrow=c(4,4))
for (i in c(3,5:16)){
  qqnorm(motif_counts_3node_5pct[,i], ylab=paste("IsoClass_",i,"_quantiles",sep=""))
  qqline(motif_counts_3node_5pct[,i])
}
dev.off()
pdf(file="~/Documents/SynCom_Modeling/param_files_enrichment/qqplots/stode/motif_countsPopulation_3node_10pct_stode.pdf")
par(mfrow=c(4,4))
for (i in c(3,5:16)){
  qqnorm(motif_counts_3node_10pct[,i], ylab=paste("IsoClass_",i,"_quantiles",sep=""))
  qqline(motif_counts_3node_10pct[,i])
}
dev.off()
pdf(file="~/Documents/SynCom_Modeling/param_files_enrichment/qqplots/stode/motif_countsPopulation_3node_15pct_stode.pdf")
par(mfrow=c(4,4))
for (i in c(3,5:16)){
  qqnorm(motif_counts_3node_15pct[,i], ylab=paste("IsoClass_",i,"_quantiles",sep=""))
  qqline(motif_counts_3node_15pct[,i])
}
dev.off()
#qqnorm plots 3node high performing PSCs counts
pdf(file="~/Documents/SynCom_Modeling/param_files_enrichment/qqplots/stode/motif_countsHP_3node_15pct_stode.pdf")
par(mfrow=c(4,4))
for (i in c(3,5:16)){
  qqnorm(motif_countsHP_3node_15pct[,i], ylab=paste("IsoClass_",i,"_quantiles",sep=""))
  qqline(motif_countsHP_3node_15pct[,i])
}
dev.off()
pdf(file="~/Documents/SynCom_Modeling/param_files_enrichment/qqplots/stode/motif_countsHP_3node_10pct_stode.pdf")
par(mfrow=c(4,4))
for (i in c(3,5:16)){
  qqnorm(motif_countsHP_3node_10pct[,i], ylab=paste("IsoClass_",i,"_quantiles",sep=""))
  qqline(motif_countsHP_3node_10pct[,i])
}
dev.off()
pdf(file="~/Documents/SynCom_Modeling/param_files_enrichment/qqplots/stode/motif_countsHP_3node_5pct_stode.pdf")
par(mfrow=c(4,4))
for (i in c(3,5:16)){
  qqnorm(motif_countsHP_3node_5pct[,i], ylab=paste("IsoClass_",i,"_quantiles",sep=""))
  qqline(motif_countsHP_3node_5pct[,i])
}
dev.off()

#below is partilly developed script for generating qqnorm plots of 4 node motif counts...havent implemented
#qqnorm plots 4node evolved PSCs counts - cant make and export these plots in bulk, too many plots.....not sure it would be useful anyway
#pdf(file="~/Documents/SynCom_Modeling/param_files_enrichment/qqplots/motif_countsEV_4node_15pct.pdf")
#par(mfrow=c(4,4))
#for (i in 4){#,8:9,13:15,17:22,25:27,30:33,36:39,41:62,64:120,122:218)){
#  qqnorm(motif_countsEV_4node_15pct[,i], ylab=paste("IsoClass_",i,"_quantiles",sep=""))
#  qqline(motif_countsEV_3node_15pct[,i])
#}
#dev.off()

write.csv(PearsonPval_motif_counts_3node,file="~/Documents/SynCom_Modeling/param_files_enrichment/stode/PearsonPval_motif_counts_3node_stode.csv")
write.csv(PearsonPval_motif_counts_4node,file="~/Documents/SynCom_Modeling/param_files_enrichment/stode/PearsonPval_motif_counts_4node_stode.csv")
#write.csv(pval_motif_counts_3node,file="~/Documents/SynCom_Modeling/param_files_enrichment/stode/pval_motif_counts_3node_enrich_stode.csv")
#write.csv(pval_motif_counts_4node,file="~/Documents/SynCom_Modeling/param_files_enrichment/stode/pval_motif_counts_4node_enrich_stode.csv")

#outputting motif counts - will use these to manually do 5 comm count window
write.csv(motif_counts_3node,file="~/Documents/SynCom_Modeling/param_files_enrichment/motif_counts/stode/motif_counts_3node.csv")
write.csv(motif_counts_3node_5pct,file="~/Documents/SynCom_Modeling/param_files_enrichment/motif_counts/stode/motif_counts_3node_5pct.csv")
write.csv(motif_counts_3node_10pct,file="~/Documents/SynCom_Modeling/param_files_enrichment/motif_counts/stode/motif_counts_3node_10pct.csv")
write.csv(motif_counts_3node_15pct,file="~/Documents/SynCom_Modeling/param_files_enrichment/motif_counts/stode/motif_counts_3node_15pct.csv")
write.csv(motif_counts_4node,file="~/Documents/SynCom_Modeling/param_files_enrichment/motif_counts/stode/motif_counts_4node.csv")
write.csv(motif_counts_4node_5pct,file="~/Documents/SynCom_Modeling/param_files_enrichment/motif_counts/stode/motif_counts_4node_5pct.csv")
write.csv(motif_counts_4node_10pct,file="~/Documents/SynCom_Modeling/param_files_enrichment/motif_counts/stode/motif_counts_4node_10pct.csv")
write.csv(motif_counts_4node_15pct,file="~/Documents/SynCom_Modeling/param_files_enrichment/motif_counts/stode/motif_counts_4node_15pct.csv")

#outputting HP community motif counts
write.csv(motif_countsHP_3node,file="~/Documents/SynCom_Modeling/param_files_enrichment/motif_counts/stode/motif_countsHP_3node.csv")
write.csv(motif_countsHP_3node_5pct,file="~/Documents/SynCom_Modeling/param_files_enrichment/motif_counts/stode/motif_countsHP_3node_5pct.csv")
write.csv(motif_countsHP_3node_10pct,file="~/Documents/SynCom_Modeling/param_files_enrichment/motif_counts/stode/motif_countsHP_3node_10pct.csv")
write.csv(motif_countsHP_3node_15pct,file="~/Documents/SynCom_Modeling/param_files_enrichment/motif_counts/stode/motif_countsHP_3node_15pct.csv")
write.csv(motif_countsHP_4node,file="~/Documents/SynCom_Modeling/param_files_enrichment/motif_counts/stode/motif_countsHP_4node.csv")
write.csv(motif_countsHP_4node_5pct,file="~/Documents/SynCom_Modeling/param_files_enrichment/motif_counts/stode/motif_countsHP_4node_5pct.csv")
write.csv(motif_countsHP_4node_10pct,file="~/Documents/SynCom_Modeling/param_files_enrichment/motif_counts/stode/motif_countsHP_4node_10pct.csv")
write.csv(motif_countsHP_4node_15pct,file="~/Documents/SynCom_Modeling/param_files_enrichment/motif_counts/stode/motif_countsHP_4node_15pct.csv")

#outputting EV community motif counts
write.csv(motif_countsEV_3node,file="~/Documents/SynCom_Modeling/param_files_enrichment/motif_counts/stode/motif_countsEV_3node.csv")
write.csv(motif_countsEV_3node_5pct,file="~/Documents/SynCom_Modeling/param_files_enrichment/motif_counts/stode/motif_countsEV_3node_5pct.csv")
write.csv(motif_countsEV_3node_10pct,file="~/Documents/SynCom_Modeling/param_files_enrichment/motif_counts/stode/motif_countsEV_3node_10pct.csv")
write.csv(motif_countsEV_3node_15pct,file="~/Documents/SynCom_Modeling/param_files_enrichment/motif_counts/stode/motif_countsEV_3node_15pct.csv")
write.csv(motif_countsEV_4node,file="~/Documents/SynCom_Modeling/param_files_enrichment/motif_counts/stode/motif_countsEV_4node.csv")
write.csv(motif_countsEV_4node_5pct,file="~/Documents/SynCom_Modeling/param_files_enrichment/motif_counts/stode/motif_countsEV_4node_5pct.csv")
write.csv(motif_countsEV_4node_10pct,file="~/Documents/SynCom_Modeling/param_files_enrichment/motif_counts/stode/motif_countsEV_4node_10pct.csv")
write.csv(motif_countsEV_4node_15pct,file="~/Documents/SynCom_Modeling/param_files_enrichment/motif_counts/stode/motif_countsEV_4node_15pct.csv")
