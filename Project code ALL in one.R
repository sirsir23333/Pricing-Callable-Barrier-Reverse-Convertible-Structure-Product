library(lubridate)
current_time1 <- Sys.time()
hour_of_day1 <- hour(current_time1)
minute_of_hour1 <- minute(current_time1)
second_of_minute1 <- second(current_time1)

## Functions
# GBM
SimGBMexact<-function(Nsim,S0,v,sigma,Deltat,T){
  m=T/Deltat # number of periods
  S=matrix(S0,nrow=Nsim,ncol=m+1)
  for(i in 1:Nsim){
    Z<-rnorm(m) # ask for m rnorm
    for(j in 2:(m+1)){ # compute value of each time step in each simulation
      S[i,j]=S[i,j-1]*exp(v*Deltat+sigma*sqrt(Deltat)*Z[j-1]) # GBM formula
    }
  }
  S
}
#AV_Sim
SimGBMexactAV<-function(Msim,S0,v,sigma,Deltat,T,collate=FALSE){
  m=T/Deltat
  S=matrix(S0,nrow=Msim,ncol=m+1)
  Stilde=matrix(S0,nrow=Msim,ncol=m+1)
  for(i in 1:Msim){
    Z<-rnorm(m)
    for(j in 2:(m+1)){
      S[i,j]=S[i,j-1]*exp(v*Deltat+sigma*sqrt(Deltat)*Z[j-1])
      Stilde[i,j]=Stilde[i,j-1]*exp(v*Deltat+sigma*sqrt(Deltat)*(-Z[j-1]))
    }
  }
  if(collate){
    out=matrix(0,2*Msim,ncol=m+1)
    for(i in 1:Msim){
      out[(2*i-1),]=S[i,]
      out[2*i,]=Stilde[i,]
    }
    return(out)
  } 
  else{return(rbind(S,Stilde))
  }
}
NotePayoff<-function(Spath,Initial=56.22,K=36.543,C=0.08,R=17.7873,D=1000){
  CouponFraction<-(440-253-i)/254
  #first coupon payment date , i=59
  if(1>=CouponFraction & CouponFraction>0.75){
    C=0.08*CouponFraction
  }
  else if(CouponFraction>0.5){
    C=0.06*(CouponFraction/0.75)
  }
  else if(CouponFraction>0.25){
    C=0.04*(CouponFraction/0.5)
  }
  else{
    C=0.02*(CouponFraction/0.25)
  }
  T=length(Spath)
  if((Spath[T]>=Initial)|(min(Spath)>K)){
    NotePay=D+C*D
  } else{
    NotePay=R*Spath[T]+C*D
  }
  #print(i)
  NotePay
}

NotePayoffNoneBarrier<-function(Spath,Initial=56.22,K=36.543,C=0.08,R=17.7873,D=1000){
  CouponFraction<-(440-253-i)/254
  if(1>=CouponFraction & CouponFraction>0.75){
    C=0.08*CouponFraction
  }
  else if(CouponFraction>0.5){
    C=0.06*(CouponFraction/0.75)
  }
  else if(CouponFraction>0.25){
    C=0.04*(CouponFraction/0.5)
  }
  else{
    C=0.02*(CouponFraction/0.25)
  }
  T=length(Spath)
  if(Spath[T]>=Initial){
    NotePay=D+C*D
  } else{
    NotePay=R*Spath[T]+C*D
  }
  
  NotePay
}

EMS<-function(SimPaths,Deltat,r){ 
  Nsim=nrow(SimPaths); m=ncol(SimPaths)-1 
  S=matrix(0,Nsim,m+1); Z=matrix(0,Nsim+1,m) 
  S[,1]=SimPaths[,1]
  for(j in 2:(m+1)){ 
    Z[1:Nsim,j-1]=S[,j-1]*SimPaths[,j]/SimPaths[,j-1] 
    Z[Nsim+1,j-1]=exp(-r*((j-1)*Deltat))*mean(Z[1:Nsim,j-1]) 
    S[,j]=SimPaths[,1]*Z[1:Nsim,j-1]/Z[Nsim+1,j-1]
  }
  S }

SimGBMpmh<-function(Nsim,S0,v,sigma,Deltat,T,h){
  m=T/Deltat # number of periods
  Splush=matrix(S0+h,nrow=Nsim,ncol=m+1)
  S=matrix(S0,nrow=Nsim,ncol=m+1)
  Sminush=matrix(S0-h,nrow=Nsim,ncol=m+1)
  for(i in 1:Nsim){
    Z<-rnorm(m)
    for(j in 2:(m+1)){
      Splush[i,j]=Splush[i,j-1]*exp(v*Deltat+sigma*sqrt(Deltat)*Z[j-1])
      S[i,j]=S[i,j-1]*exp(v*Deltat+sigma*sqrt(Deltat)*Z[j-1])
      Sminush[i,j]=Sminush[i,j-1]*exp(v*Deltat+sigma*sqrt(Deltat)*Z[j-1])
    }
  }
  list(Splush=Splush,S=S,Sminush=Sminush)
}

## Project code
# Load hist data and compute anuual mean and sd
LOGN<-read.csv("LOGN.SW_non.csv",header=TRUE) # read data 11 Aug 2022 - 9 Nov 2023
rates_termstructure<-read.csv("Rates with term structure.csv",header=TRUE)
dt=1/252
LOGNprices<-as.numeric(as.vector(LOGN$Adj.Close))
LOGNprices<-LOGNprices[!is.na(LOGNprices)]

# Parameters
Nsim <- 10000
Msim <- 5000
T <- 188  # number of trading days in UK between 10/08/23 and 31/07/24
r2<- read.csv("Rates.csv",header=TRUE)


# Pricing
results <- matrix(nrow = 66, ncol = 5)
colnames(results) <- c("i", "NoteValue", "NoteValueVar", "NoteValueLow", "NoteValueUp")
resultsAV <- matrix(nrow = 66, ncol = 5)
colnames(resultsAV) <- c("i", "NoteValueAV", "NoteValueVarAV", "NoteValueLowAV", "NoteValueUpAV")
resultsCV <- matrix(nrow = 66, ncol = 5)
colnames(resultsCV) <- c("i", "NoteValueCV", "NoteValueVaCVr", "NoteValueLowCV", "NoteValueUpCV")
resultsEMS <-matrix(nrow=66,ncol=5)
colnames(resultsEMS) <- c("i", "NoteValueCV", "NoteValueVaCVr", "NoteValueLowCV", "NoteValueUpCV")
resultsCall<-matrix(nrow=66,ncol=5)
colnames(resultsCall) <- c("i", "NoteValueCall", "NoteValueVarCall", "NoteValueLowCall", "NoteValueUpCall")
resultsCallCV<-matrix(nrow=66,ncol=5)
colnames(resultsCallCV) <- c("i", "NoteValueCallCV", "NoteValueVarCallCV", "NoteValueLowCallCV", "NoteValueUpCallCV")
results_termR <- matrix(nrow = 66, ncol = 5)
colnames(results_termR) <- c("i", "NoteValue_termR", "NoteValueVar_termR", "NoteValueLow_termR", "NoteValueUp_termR")
resultsNoneBarrier<-matrix(nrow=66,ncol=5)
colnames(resultsNoneBarrier) <- c("i", "NoteValueNonBarrier", "NoteValueVarNonBarrier", "NoteValueLowNonBarrier", "NoteValueUpNonBarrier")

resultsNoCallDelta <- matrix(nrow = 66, ncol = 2)
colnames(resultsNoCallDelta) <- c("i", "Delta")
resultsNoCallGamma <- matrix(nrow = 66, ncol = 2)
colnames(resultsNoCallGamma) <- c("i", "Gamma")
resultsCallDelta <- matrix(nrow = 66, ncol = 2)
colnames(resultsCallDelta) <- c("i", "Delta")
resultsCallGamma <- matrix(nrow = 66, ncol = 2)
colnames(resultsCallGamma) <- c("i", "Gamma")

for (i in 0:65) {
  if(60-i<=0){
    Exdates<-c(120-i)
  }else{
    Exdates<-c(60-i,120-i)
  }
  Extimes<-length(Exdates)
  ExtimesForGreeks<-length(Exdates)
  Exdates<-Exdates[order(Exdates)] # sort exercise dates indices in ascending order Extimes=length(Exdates)
  if(Exdates[Extimes]==T/dt){Exdates=Exdates[-Extimes]} # terminal date should be excluded
  N1<-100
  histdata <- LOGNprices[1 + i:252 + i]  # select the historical data for pricing
  discount_factor <- rates_termstructure[(1+i),4]
  r <- r2[(1+i),8]*0.01
  n0 <- length(histdata)
  logprices <- log(histdata)
  logreturns <- logprices[2:n0] - logprices[1:(n0 - 1)]  # ln(P_t) - ln(P_t-1)
  v <- mean(logreturns) / dt  # annual mean
  sigma <- sd(logreturns) / sqrt(dt)  # annual sd
  
  # risk-neutral valuation (standard MC)
  Tminust <- (T - i) / 252
  St=histdata[n0]
  v=r-(sigma^2)/2
  set.seed(4518)
  SimLOGN <- SimGBMexact(Nsim, St, v, sigma, dt, Tminust)
  price_before_t<-LOGNprices[184:252+i]
  #price_before_t
  
  Histdata<-matrix(rep(price_before_t,Nsim),ncol=length(price_before_t),byrow=T)
  wholedata<-cbind(Histdata[,-ncol(Histdata)],SimLOGN)
  SimNotes <- exp(-r * Tminust) * apply(wholedata, 1, NotePayoff)
  #print(SimNotes)
  NoteValue <- mean(SimNotes)
  #print(NoteValue)
  NoteValueVar <- var(SimNotes) / Nsim
  NoteValueLow <- NoteValue - 1.96 * sqrt(NoteValueVar)
  NoteValueUp <- NoteValue + 1.96 * sqrt(NoteValueVar)
  #print(NoteValue)
  
  SimLOGNAV <- SimGBMexactAV(Msim, St, v, sigma, dt, Tminust)
  HistdataAV<-matrix(rep(price_before_t,2*Msim),ncol=length(price_before_t),byrow=T)
  wholedataAV<-cbind(HistdataAV[,-ncol(HistdataAV)],SimLOGNAV)
  SimNotesAV <- exp(-r * Tminust) * apply(wholedataAV, 1, NotePayoff)
  NoteValueAV <- mean(SimNotesAV)
  NoteValueVarAV <- var(SimNotesAV) / Nsim
  NoteValueLowAV <- NoteValueAV - 1.96 * sqrt(NoteValueVarAV)
  NoteValueUpAV <- NoteValueAV + 1.96 * sqrt(NoteValueVarAV)
  
  SimLOGNCV <- SimGBMexact(Nsim, St, v, sigma, dt, Tminust)
  Histdata<-matrix(rep(price_before_t,Nsim),ncol=length(price_before_t),byrow=T)
  Histdatapre<-matrix(rep(price_before_t,N1),ncol=length(price_before_t),byrow=T)
  wholedata<-cbind(Histdata[,-ncol(Histdata)],SimLOGNCV)
  wholedatapre<-cbind(Histdatapre[,-ncol(Histdatapre)],SimLOGNCVpre)
  SimNotespre<-exp(-r*Tminust)*apply(wholedatapre,1,NotePayoff)
  c=-cov(SimNotespre,wholedatapre[,ncol(wholedatapre)])/(St^2*exp(2*r*Tminust)*(exp(sigma^2*Tminust)-1))
  SimNotesCV<-exp(-r*Tminust)*apply(wholedata,1,NotePayoff)+c*(SimLOGNCV[,ncol(SimLOGNCV)]-St*exp(r*Tminust))
  NoteCVValue=mean(SimNotesCV)
  NoteCVValueVar=var(SimNotesCV)/Nsim 
  
  SimLOGN <- SimGBMexact(Nsim, St, v, sigma, dt, Tminust)
  EMSpaths<-EMS(SimLOGN,dt,r)
  EMSpaths<-EMSpaths[,-1]
  Histdata<-matrix(rep(price_before_t,Nsim),ncol=length(price_before_t),byrow=T)
  wholedataEMS<-cbind(Histdata[,-ncol(Histdata)],EMSpaths)
  SimNotesEMS <- exp(-r * Tminust) * apply(wholedataEMS, 1, NotePayoff)
  NoteValueEMS <- mean(SimNotesEMS)
  NoteValueVarEMS <- var(SimNotesEMS) / Nsim
  NoteValueLowEMS <- NoteValueEMS - 1.96 * sqrt(NoteValueVarEMS)
  NoteValueUpEMS <- NoteValueEMS + 1.96 * sqrt(NoteValueVarEMS)
  
  SimLOGN_termR <- SimGBMexact(Nsim, St, v, sigma, dt, Tminust)
  price_before_t<-LOGNprices[184:(252+i)]
  Histdata_termR<-matrix(rep(price_before_t,Nsim),ncol=length(price_before_t),byrow=T)
  wholepath_termR<-cbind(Histdata_termR[,-ncol(Histdata_termR)],SimLOGN_termR)
  SimNotes_termR <- discount_factor * apply(wholepath_termR, 1, NotePayoff)
  NoteValue_termR<- mean(SimNotes_termR)
  NoteValueVar_termR <- var(SimNotes_termR) / Nsim
  NoteValueLow_termR <- NoteValue_termR - 1.96 * sqrt(NoteValueVar_termR)
  NoteValueUp_termR <- NoteValue_termR + 1.96 * sqrt(NoteValueVar_termR)
  
  SimLOGNNonBarrier<-SimGBMexact(Nsim, St, v, sigma, dt, Tminust)
  HistdataNonBarrier<-matrix(rep(price_before_t,Nsim),ncol=length(price_before_t),byrow=T)
  wholedataNonBarrier<-cbind(HistdataNonBarrier[,-ncol(HistdataNonBarrier)],SimLOGNNonBarrier)
  SimNotesNonBarrier <- exp(-r * Tminust) * apply(wholedataNonBarrier, 1, NotePayoffNoneBarrier)
  NoteValueNonBarrier <- mean(SimNotesNonBarrier)
  NoteValueVarNonBarrier <- var(SimNotesNonBarrier) / Nsim
  NoteValueLowNonBarrier <- NoteValueNonBarrier - 1.96 * sqrt(NoteValueVarNonBarrier)
  NoteValueUpNonBarrier <- NoteValueNonBarrier + 1.96 * sqrt(NoteValueVarNonBarrier)
  # Set value of h for Delta and Gamma
  h=St * 0.001
  
  # Simulate Delta and Gamma without callable feature
  set.seed(4518)
  SimPlusMinushLOGN<-SimGBMpmh(Nsim,St,v,sigma,dt,Tminust,h)
  SimPlusMinushLOGNSplush<-SimPlusMinushLOGN$Splush
  SimPlusMinushLOGNS<-SimPlusMinushLOGN$S
  SimPlusMinushLOGNSminush<-SimPlusMinushLOGN$Sminush
  
  SimNotesSplush<-apply(SimPlusMinushLOGNSplush,1,NotePayoff)
  SimNotesS<-apply(SimPlusMinushLOGNS,1,NotePayoff)
  SimNotesSminush<-apply(SimPlusMinushLOGNSminush,1,NotePayoff)
  NoteValueSplush=exp(-r*1)*mean(SimNotesSplush)
  NoteValueS=exp(-r*1)*mean(SimNotesS)
  NoteValueSminush=exp(-r*1)*mean(SimNotesSminush)
  
  NoteDelta=(NoteValueSplush-NoteValueSminush)/(2*h)
  NoteGamma=(NoteValueSplush-2*NoteValueS+NoteValueSminush)/(h^2)
  
  # Simulate Call
  # Step 2
  SimLOGNCall <- SimGBMexact(Nsim, St, v, sigma, dt, Tminust)
  HistdataCall<-matrix(rep(price_before_t,Nsim),ncol=length(price_before_t),byrow=T)
  wholedataCall<-cbind(HistdataCall[,-ncol(HistdataCall)],SimLOGNCall)
  m=ncol(wholedataCall)
  Extimes<-Extimes+1
  V=matrix(0,nrow=Nsim,ncol=Extimes)
  V[,Extimes]=apply(as.matrix(wholedataCall[,m]),1,NotePayoff)
  
  #Call CV
  SimLOGNCallCV <- SimGBMexact(Nsim, St, v, sigma, dt, Tminust)
  HistdataCallCV<-matrix(rep(price_before_t,Nsim),ncol=length(price_before_t),byrow=T)
  wholedataCallCV<-cbind(HistdataCallCV[,-ncol(HistdataCallCV)],SimLOGNCallCV)
  
  SimLOGNCallCVpre<-SimGBMexact(N1, St, v, sigma, dt, Tminust)
  HistdataCallCVpre<-matrix(rep(price_before_t,N1),ncol=length(price_before_t),byrow=T)
  wholedataCallCVpre<-cbind(HistdataCallCVpre[,-ncol(HistdataCallCVpre)],SimLOGNCallCVpre)
  
  V_CV=matrix(0,nrow=Nsim,ncol=Extimes)
  V_CVpre=matrix(0,nrow=100,ncol=Extimes)
  
  V_CV[,Extimes]=apply(as.matrix(wholedataCallCV[,m]),1,NotePayoff)
  V_CVpre[,Extimes]=apply(as.matrix(wholedataCallCVpre[,m]),1,NotePayoff)
  
  previousstep<-m
  for(j in Exdates[order(Exdates,decreasing=TRUE)]){
    currentstep<-j
    # Step 3a
    Extimes<-Extimes-1
    
    Y=exp(-r*(previousstep-currentstep)*dt)*V[,Extimes+1];
    X=wholedataCall[,currentstep]
    X2=X^2
    Reg<-lm(Y~X+X2)
    
    Y_CVpre=exp(-r*(previousstep-currentstep)*dt)*V_CVpre[,Extimes+1]
    X_CVpre=wholedataCallCVpre[,currentstep]
    X_CV2pre=X_CVpre^2
    #print(c(length(X_CVpre)))
    RegCVpre<-lm(Y_CVpre~X_CVpre+X_CV2pre)

    
    Y_CV=exp(-r*(previousstep-currentstep)*dt)*V_CV[,Extimes+1]
    X_CV=wholedataCallCV[,currentstep]
    X_CV2=X_CV^2
    RegCV<-lm(Y_CV~X_CV+X_CV2)
    #RegResults[,Extimes]=Reg$coefficients
    # Step 3b
    if(Extimes==2){
      redemptioncost<-1000*(1+0.08*(3/4))
    }else if(Extimes==1){
      redemptioncost<-1000*(1+0.08*(1/2))
    }else{
      redemptioncost<-1000*(1+0.08*(1/4))
    }
    
    ExerciseCond<-(redemptioncost<predict(Reg,as.data.frame(wholedataCall[,currentstep])))
    Exerciseindex<-which(ExerciseCond);
    Redemptionindex<-which(!ExerciseCond)
    V[Exerciseindex,Extimes]=redemptioncost
    
    ExerciseCondCVpre<-(redemptioncost<predict(RegCVpre,as.data.frame(wholedataCallCVpre[,currentstep])))
    ExerciseindexCVpre<-which(ExerciseCondCVpre);
    RedemptionindexCVpre<-which(!ExerciseCondCVpre)
    V_CVpre[ExerciseindexCVpre,Extimes]=redemptioncost
    
    ExerciseCond_CV<-(redemptioncost<predict(RegCV,as.data.frame(wholedataCallCV[,currentstep])))
    Exerciseindex_CV<-which(ExerciseCond_CV);
    Redemptionindex_CV<-which(!ExerciseCond_CV)
    V_CV[Exerciseindex_CV,Extimes]=redemptioncost
    
    V[Redemptionindex,Extimes]=exp(-r*(previousstep-currentstep)*dt)*V[Redemptionindex,Extimes+1]
    V_CV[Redemptionindex_CV,Extimes]=exp(-r*(previousstep-currentstep)*dt)*V[Redemptionindex_CV,Extimes+1]
    V_CVpre[RedemptionindexCVpre,Extimes]=exp(-r*(previousstep-currentstep)*dt)*V_CVpre[RedemptionindexCVpre,Extimes+1]
    
    previousstep<-j
  }
  NoteValueCall=exp(-r*(previousstep-0)*dt)*mean(V[,1]);
  SimNotesCall=exp(-r*(previousstep-0)*dt)*V[,1]
  NoteValueVarCall <- var(SimNotesCall) / Nsim
  NoteValueLowCall <- NoteValueCall - 1.96 * sqrt(NoteValueVarCall)
  NoteValueUpCall <- NoteValueCall + 1.96 * sqrt(NoteValueVarCall)
  resultsCall[i+1,] <- c(i, NoteValueCall, NoteValueVarCall, NoteValueLowCall, NoteValueUpCall)
  
  NoteValueCallCVpre=exp(-r*(previousstep-0)*dt)*mean(V_CVpre[,1]);
  SimNotesCallCVpre=exp(-r*(previousstep-0)*dt)*V_CVpre[,1]
  c_Call=-cov(SimNotesCallCVpre,wholedataCallCVpre[,ncol(wholedataCallCVpre)])/(St^2*exp(2*r*(previousstep))*(exp(sigma^2*previousstep)-1))
  
  NoteValueCallCV<-exp(-r*(previousstep-0)*dt)*mean(V_CV[,1])+c_Call*(mean(SimLOGNCallCV[,ncol(SimLOGNCallCV)]))-St*exp(r*(previousstep-0)*dt)
  
  SimNotesCallCV=exp(-r*(previousstep-0)*dt)*V_CV[,1]+c_Call*(SimLOGNCallCV[,ncol(SimLOGNCallCV)]-St*exp(r*(previousstep-0)*dt))
  
  NoteValueVarCallCV <- var(SimNotesCallCV) / Nsim
  #print(NoteValueVarCallCV)
  # Try Simulating Delta and Gamma with Callable Features
  
  # Simulate Call for S, S + h, S - h
  set.seed(4518)
  SimPlusMinushLOGNCall <- SimGBMpmh(Nsim,St,v,sigma,dt,Tminust,h)
  
  SimPlusMinushLOGNCallSplush<-SimPlusMinushLOGNCall$Splush
  SimPlusMinushLOGNCallS<-SimPlusMinushLOGNCall$S
  SimPlusMinushLOGNCallSminush<-SimPlusMinushLOGNCall$Sminush
  mSplush<-ncol(SimPlusMinushLOGNCallSplush)
  mS<-ncol(SimPlusMinushLOGNCallS)
  mSminush<-ncol(SimPlusMinushLOGNCallSminush)
  
  # Finding payoff for all 3 price path
  ExtimesForGreeks<-ExtimesForGreeks+1
  V_S=matrix(0,nrow=Nsim,ncol=ExtimesForGreeks); 
  V_Splush=matrix(0,nrow=Nsim,ncol=ExtimesForGreeks);
  V_Sminush=matrix(0,nrow=Nsim,ncol=ExtimesForGreeks);
  
  #RegResults=matrix(0,nrow=3,ncol=ExtimesForGreeks-1)
  # Step 2
  V_S[,ExtimesForGreeks]=apply(as.matrix(SimPlusMinushLOGNCallS[,mS]),1,NotePayoff) 
  V_Splush[,ExtimesForGreeks]=apply(as.matrix(SimPlusMinushLOGNCallSplush[,mSplush]),1,NotePayoff)
  V_Sminush[,ExtimesForGreeks]=apply(as.matrix(SimPlusMinushLOGNCallSminush[,mSminush]),1,NotePayoff)
  
  previousstep<-mS
  for(j in Exdates[order(Exdates,decreasing=TRUE)]){
    currentstep<-j
    # Step 3a
    ExtimesForGreeks<-ExtimesForGreeks-1
    Y_S=exp(-r*(previousstep-currentstep)*dt)*V_S[,ExtimesForGreeks+1]; 
    Y_Sminush=exp(-r*(previousstep-currentstep)*dt)*V_Sminush[,ExtimesForGreeks+1];
    Y_Splush=exp(-r*(previousstep-currentstep)*dt)*V_Splush[,ExtimesForGreeks+1];
    
    # For S
    X_S=SimPlusMinushLOGNCallS[,currentstep]
    X2_S= (X_S)^2
    Reg_S<-lm(Y_S~X_S+X2_S)
    
    # For S + h
    X_Splush=SimPlusMinushLOGNCallSplush[,currentstep]
    X2_Splush= (X_Splush)^2
    Reg_Splush<-lm(Y_Splush~X_Splush+X2_Splush)
    
    # For S - h
    X_Sminush=SimPlusMinushLOGNCallSminush[,currentstep]
    X2_Sminush= (X_Sminush)^2
    Reg_Sminush<-lm(Y_Sminush~X_Sminush+X2_Sminush)
    
    #RegResults[,ExtimesForGreeks]=Reg$coefficients
    # Step 3b
    if(Extimes==2){
      redemptioncost<-1000*(1+0.08*(3/4))
    }else if(Extimes==1){
      redemptioncost<-1000*(1+0.08*(1/2))
    }else{
      redemptioncost<-1000*(1+0.08*(1/4))
    }
    
    # Need to do for S, S+h and S-h
    #print(c(dim(SimPlusMinushLOGNCallS),currentstep))
    ExerciseCondS<-(redemptioncost<predict(Reg_S,as.data.frame(SimPlusMinushLOGNCallS[,currentstep])))
    ExerciseCondSminush<-(redemptioncost<predict(Reg_Sminush,as.data.frame(SimPlusMinushLOGNCallSplush[,currentstep])))
    ExerciseCondSplush<-(redemptioncost<predict(Reg_Splush,as.data.frame(SimPlusMinushLOGNCallSminush[,currentstep])))
    
    ExerciseindexS<-which(ExerciseCondS); 
    ExerciseindexSminush<-which(ExerciseCondSminush); 
    ExerciseindexSplush<-which(ExerciseCondSplush); 
    
    RedemptionindexS<-which(!ExerciseCondS)
    RedemptionindexSminush<-which(!ExerciseCondSminush)
    RedemptionindexSplush<-which(!ExerciseCondSplush)
    
    V_S[ExerciseindexS,ExtimesForGreeks]=redemptioncost
    V_Sminush[ExerciseindexSminush,ExtimesForGreeks]=redemptioncost
    V_Splush[ExerciseindexSplush,ExtimesForGreeks]=redemptioncost
    
    V_S[RedemptionindexS,ExtimesForGreeks]=exp(-r*(previousstep-currentstep)*dt)*V_S[RedemptionindexS,ExtimesForGreeks+1]
    V_Sminush[RedemptionindexSminush,ExtimesForGreeks]=exp(-r*(previousstep-currentstep)*dt)*V_Sminush[RedemptionindexSminush,ExtimesForGreeks+1]
    V_Splush[RedemptionindexSplush,ExtimesForGreeks]=exp(-r*(previousstep-currentstep)*dt)*V_Splush[RedemptionindexSplush,ExtimesForGreeks+1]
    
    previousstep<-j
  }
  
  SimNotesCallS=exp(-r*(previousstep-0)*dt)*V_S[,1]
  SimNotesCallSminush=exp(-r*(previousstep-0)*dt)*V_Sminush[,1]
  SimNotesCallSplush=exp(-r*(previousstep-0)*dt)*V_Splush[,1]
  
  NoteValueCallS <- mean(SimNotesCallS)
  NoteValueCallSminush <- mean(SimNotesCallSminush)
  NoteValueCallSplush <- mean(SimNotesCallSplush)
  
  NoteDeltaCall=(NoteValueCallSplush-NoteValueCallSminush)/(2*h)
  NoteGammaCall=(NoteValueCallSplush-2*NoteValueCallS+NoteValueCallSminush)/(h^2)
  
  #plot(RegResults[1,],type="l",ylim=c(min(RegResults),max(RegResults))) 
  #lines(RegResults[2,],col="red")
  # store results
  resultsAV[i + 1, ] <- c(i, NoteValueAV, NoteValueVarAV, NoteValueLowAV, NoteValueUpAV)
  resultsCV[i + 1, ] <- c(i, NoteCVValue, NoteCVValueVar, NoteCVValue-1.96*sqrt(NoteCVValueVar), NoteCVValue+1.96*sqrt(NoteCVValueVar))
  #print(c(NoteValueCallCV,NoteValueVarCallCV))
  
  resultsCallCV[i+1,]<-c(i,NoteValueCallCV,NoteValueVarCallCV,NoteValueCallCV-1.96*sqrt(NoteValueVarCallCV),NoteValueCallCV+1.96*sqrt(NoteValueVarCallCV))

  results[i + 1, ] <- c(i, NoteValue, NoteValueVar, NoteValueLow, NoteValueUp)
  resultsEMS[i+1,]<-c(i,NoteValueEMS,NoteValueVarEMS,NoteValueLowEMS,NoteValueUpEMS)
  resultsCall[i+1,]<-c(i,NoteValueCall,NoteValueVarCall,NoteValueLowCall,NoteValueUpCall)
  results_termR[i + 1, ] <- c(i, NoteValue_termR, NoteValueVar_termR, NoteValueLow_termR, NoteValueUp_termR)
  resultsNoneBarrier[i+1,]<-c(i, NoteValueNonBarrier, NoteValueVarNonBarrier, NoteValueLowNonBarrier, NoteValueUpNonBarrier)
  resultsNoCallDelta[i + 1, ] <- c(i, NoteDelta)
  resultsNoCallGamma[i + 1, ] <- c(i, NoteGamma)
  resultsCallDelta[i + 1, ] <- c(i, NoteDeltaCall)
  resultsCallGamma[i + 1, ] <- c(i, NoteGammaCall)
}

# Visualisation (model vs actual)
# results
actual<-read.csv("non_call_option.csv",header=TRUE)
actual_price= actual$Actual.price[1:66]
ModelandActual<-cbind(resultsCV,actual[-nrow(actual),])
actualdata <- ModelandActual[, 7]
modeldataCV <- ModelandActual[, 2]
modeldataMC <- results[,2]
modeldataAV <- resultsAV[,2]
modeldataEMS<-resultsEMS[,2]
modeldataCall<-resultsCall[,2]
modeldataCallCV<-resultsCallCV[,2]
modeldata_termR<-results_termR[,2]
modeldataNonBarrier<-resultsNoneBarrier[,2]

#These are all the lines for our report, we can have different combinations.
plot(actualdata, col = "blue", type = "l", lty = 1, ylim = range(c(900,1035)), ylab = "Structure price", xlab = "Number of trading days from 09/08/23", main = "Model price VS Actual price")
lines(modeldataCV, col = "red", type = "l", lty = 1)
lines(modeldataMC, col = "green", type = "l", lty = 1)
lines(modeldataAV, col = "black", type = "l", lty = 1)
lines(modeldataEMS,col="yellow",type='l',lty=1)
lines(modeldataCall,col='purple',type='l',lty=1)
lines(modeldata_termR,col='grey',type='l',lty=1)
lines(modeldataNonBarrier,col='orange',type='l',lty=1)
lines(modeldataCallCV,col='aliceblue',type='l',lty=1)

#Add any legends you want, below are two examples
legend("bottomright", legend = c("Actual Non-Call Price", "Model CV","Model Non-Callable","Model AV","Model EMS",'Model Call', "Model Non-Barrier"), col = c("blue", "red","green","black","yellow",'purple','orange'), lty = 1:1, cex = 0.8)  
legend("bottomright", legend = c("Actual Non-Call Price", "Model Barrier", "Model Non-Barrier"), col = c("blue", "red","green","black","yellow",'purple','orange'), lty = 1:1, cex = 0.8)  

# Plot Delta and Gamma without Call
modeldataNoCallDelta <- resultsNoCallDelta[, 2]
modeldataNoCallGamma <- resultsNoCallGamma[, 2]
plot(modeldataNoCallDelta, col = "purple", type = "l", lty = 1, ylim = range(modeldataNoCallDelta), ylab = "Delta", xlab = "Number of trading days from 09/08/23", main = "Delta Without Callable Feature VS Trading Days")
plot(modeldataNoCallGamma, col = "brown", type = "l", lty = 1, ylim = range(modeldataNoCallGamma), ylab = "Gamma", xlab = "Number of trading days from 09/08/23", main = "Gamma Without Callable Feature VS Trading Days")

# Plot Delta and Gamma with Call
modeldataCallDelta <- resultsCallDelta[, 2]
modeldataCallGamma <- resultsCallGamma[, 2]
plot(modeldataCallDelta, col = "pink", type = "l", lty = 1, ylim = range(modeldataCallDelta), ylab = "Delta", xlab = "Number of trading days from 09/08/23", main = "Delta With Callable Feature VS Trading Days")
plot(modeldataCallGamma, col = "grey", type = "l", lty = 1, ylim = range(modeldataCallGamma), ylab = "Gamma", xlab = "Number of trading days from 09/08/23", main = "Gamma With Callable Feature VS Trading Days")

# Can try to plot the Delta and Gamma with and without call together
plot(modeldataNoCallDelta, col = "purple", type = "l", lty = 1, ylim = range(c(modeldataNoCallDelta,modeldataCallDelta)), ylab = "Delta", xlab = "Number of trading days from 09/08/23", main = "Callable vs NonCallable For Delta")
lines(modeldataCallDelta, col = "pink", type = "l", lty = 1)
legend("topright", legend = c("Non Callable", "Callable"), col = c("purple", "pink"), lty = 1:1, cex = 0.8)  

plot(modeldataNoCallGamma, col = "brown", type = "l", lty = 1, ylim = range(c(modeldataNoCallGamma,modeldataCallGamma)), ylab = "Gamma", xlab = "Number of trading days from 09/08/23", main = "Callable vs NonCallable For Gamma")
lines(modeldataCallGamma, col = "grey", type = "l", lty = 1)
legend("topright", legend = c("Non Callable", "Callable"), col = c("brown", "grey"), lty = 1:1, cex = 0.8)

print(mean(modeldataCall))
print(mean(modeldataMC))
## mean of variance 
MeanVarMC=mean(results[,3])
MeanVarAV=mean(resultsAV[,3])
MeanVarCV=mean(resultsCV[,3])
MeanVarEMS=mean(resultsEMS[,3])
MeanVarCall=mean(resultsCall[,3])
MeanVarNonBarrier=mean(resultsNoneBarrier[,3])
MeanVarCallCV=mean(resultsCallCV[,3])
print(MeanVarMC)
print(MeanVarAV)
print(MeanVarCV)
print(MeanVarEMS)
print(MeanVarCall)
print(MeanVarCallCV)

current_time2 <- Sys.time()
hour_of_day2 <- hour(current_time2)
minute_of_hour2 <- minute(current_time2)
second_of_minute2 <- second(current_time2)
print(paste("Current time1:", hour_of_day1, "hours", minute_of_hour1, "minutes", second_of_minute1, "seconds"))

print(paste("Current time2:", hour_of_day2, "hours", minute_of_hour2, "minutes", second_of_minute2, "seconds"))


