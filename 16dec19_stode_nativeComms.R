library(rootSolve)
n<-1
num_of_comms <-250
final_parameters <- matrix(ncol=30, nrow=num_of_comms)
final_states <- matrix(ncol=5,nrow=num_of_comms)
final_eigs <- matrix(ncol=5,nrow=num_of_comms)
# distribution for interaction coefficients  
param_mean <- (-0.375)
param_sd <- 1
# distribution for mu (growth rates)
mu_mean <- 0.002
mu_sd <- 0.0002
# assigns initial conditions, in this case all species start at same abundance, changes as the comm evolves
start_amt<-0.1

#begin native Community generation

for (n in 1:num_of_comms){
repeat{
#using solve() to solve linearized gLV - solutions from which will be the initial guess provided to nonlinear solver 'stode' 
#if statement within the repeat loop which starts below (the one using solve() to find roots of linearized gLV), ensures only solutions containing positive numbers are passed to stode
  repeat{
    
    muX1<-abs(rnorm(1,mu_mean,sd=mu_sd))
    muX2<-abs(rnorm(1,mu_mean,sd=mu_sd)) 
    muX3<-abs(rnorm(1,mu_mean,sd=mu_sd))
    muX4<-abs(rnorm(1,mu_mean,sd=mu_sd))
    muX5<-abs(rnorm(1,mu_mean,sd=mu_sd))
    a11<-abs(rnorm(1,mean=param_mean,sd=param_sd))*-1
    a12<-rnorm(1,mean=param_mean,sd=param_sd)
    a13<-rnorm(1,param_mean,sd=param_sd)
    a14<-rnorm(1,param_mean,sd=param_sd)
    a15<-rnorm(1,param_mean,sd=param_sd)
    a21<-rnorm(1,param_mean,sd=param_sd)
    a22<-abs(rnorm(1,param_mean,sd=param_sd))*-1
    a23<-rnorm(1,param_mean,sd=param_sd)
    a24<-rnorm(1,param_mean,sd=param_sd)
    a25<-rnorm(1,param_mean,sd=param_sd)
    a31<-rnorm(1,param_mean,sd=param_sd)
    a32<-rnorm(1,param_mean,sd=param_sd)
    a33<-abs(rnorm(1,param_mean,sd=param_sd))*-1
    a34<-rnorm(1,param_mean,sd=param_sd)
    a35<-rnorm(1,param_mean,sd=param_sd)
    a41<-rnorm(1,param_mean,sd=param_sd)
    a42<-rnorm(1,param_mean,sd=param_sd)
    a43<-rnorm(1,param_mean,sd=param_sd)
    a44<-abs(rnorm(1,param_mean,sd=param_sd))*-1
    a45<-rnorm(1,param_mean,sd=param_sd)
    a51<-rnorm(1,param_mean,sd=param_sd)
    a52<-rnorm(1,param_mean,sd=param_sd)
    a53<-rnorm(1,param_mean,sd=param_sd)
    a54<-rnorm(1,param_mean,sd=param_sd)
    a55<-abs(rnorm(1,param_mean,sd=param_sd))*-1
    
    At=matrix(c(a11,a12,a13,a14,a15,a21,a22,a23,a24,a25,a31,a32,a33,a34,a35,a41,a42,a43,a44,a45,a51,a52,a53,a54,a55), nrow=5, ncol=5)
    A=t(At)
    b=matrix(c(-muX1,-muX2,-muX3,-muX4,-muX5),nrow=5,ncol=1)
    x<-solve(A,b)
    if(x[1]>0&x[2]>0&x[3]>0&x[4]>0&x[5]>0){
    zerocheck<-A%*%x-b
    #print(zerocheck)
    break
    }
  }

  #linearized 'guess' of solution from solve() being supplied to nonlinear solver as initial guess
  state<-c(X1=x[1],X2=x[2],X3=x[3],X4=x[4],X5=x[5])
  
# assigns growth rates for all species and interaction coefficients as parameters for nonlinear solver
  parametersNATCOM<-c(muX1,
                      muX2, 
                      muX3,
                      muX4,
                      muX5,
                      a11,
                      a12,
                      a13,
                      a14,
                      a15,
                      a21,
                      a22,
                      a23,
                      a24,
                      a25,
                      a31,
                      a32,
                      a33,
                      a34,
                      a35,
                      a41,
                      a42,
                      a43,
                      a44,
                      a45,
                      a51,
                      a52,
                      a53,
                      a54,
                      a55)
  # function for generating differential equations according to the Lotka-Volterra model
  NonLinearGLV<-function(t, state, parameters) {
    with(as.list(c(state, parameters)),{
      dX1 <- muX1*X1 + a11*X1*X1 + a12*X1*X2 + a13*X1*X3 + a14*X1*X4 + a15*X1*X5
      dX2 <- muX2*X2 + a22*X2*X2 + a21*X2*X1 + a23*X2*X3 + a24*X2*X4 + a25*X2*X5
      dX3 <- muX3*X3 + a33*X3*X3 + a31*X3*X1 + a32*X3*X2 + a34*X3*X4 + a35*X3*X5
      dX4 <- muX4*X4 + a44*X4*X4 + a41*X4*X1 + a42*X4*X2 + a43*X4*X3 + a45*X4*X5
      dX5 <- muX5*X5 + a55*X5*X5 + a51*X5*X1 + a52*X5*X2 + a53*X5*X3 + a54*X5*X4
      list(c(dX1,dX2,dX3,dX4,dX5))
    })
  }
  out_4<-stode(y=state,fun=NonLinearGLV,parms=parametersNATCOM, pos=TRUE)
  jac<-jacobian.full(y=c(out_4$y[1],out_4$y[2],out_4$y[3],out_4$y[4],out_4$y[5]),func=NonLinearGLV, parms=parametersNATCOM)
  eig<-eigen(jac)
#ensures eigenvalues are real and negative. If so, prints final states and parameters  
if(!is.complex(eig$values[1])&!is.complex(eig$values[2])&!is.complex(eig$values[3])&!is.complex(eig$values[4])&!is.complex(eig$values[5])){
if(eig$values[1]<0&eig$values[2]<0&eig$values[3]<0&eig$values[4]<0&eig$values[5]<0){
  final_states[n,]=out_4$y
  final_parameters[n,]<-parametersNATCOM
  final_eigs[n,]<-eig$values
break    
  }
}
  
}
n=n+1
print(n)
}
write.table(final_states, file = "~/Documents/SynCom_Modeling/NativeCommunities_stode_final_states.csv", sep = ",", col.names = F, row.names = F, append = TRUE)
write.table(final_parameters, file = "~/Documents/SynCom_Modeling/NativeCommunities_stode_final_parameters.csv", sep = ",", col.names = F, row.names = F, append = TRUE)
write.table(final_eigs, file = "~/Documents/SynCom_Modeling/NativeCommunities_stode_final_eigs.csv", sep = ",", col.names = F, row.names = F, append = TRUE)