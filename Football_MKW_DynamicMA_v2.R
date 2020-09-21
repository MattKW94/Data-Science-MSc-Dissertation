library(skellam)
library(lattice)

# setwd("C:/Users/matth/Dropbox/Matthew Dissertation/R Work/Football_Data")

# Ranked probability score function
RPS=function(X,P){
  PP=t(apply(P,1,cumsum))[,1:2]
  XX=t(apply(X,1,cumsum))[,1:2]
  RES=apply((PP-XX)^2,1,mean)
  return(RES)}

# Data needs fixing for 1996 and 2020 seasons
I=diag(rep(1,3)) # I matrix, for calculating metrics
ODO=rep(0,25) # Initial ODOs, ?
SEA=as.character(2006:2020) # Seasons years
NS=paste(SEA,".csv",sep="") # Names of seasons files
PICS=paste(SEA,".png",sep="") # Names of seasons pic files
PICSpdf=paste(SEA,".pdf",sep="") # Names of seasons pic pdf files
All=as.list(1:25) # Attacking and defensive strengths
PRES=as.list(1:25) # Pearson residuals
DRES=as.list(1:25) # Deviance residuals
NAMESIN=as.list(1:25) # NAMES in each season
DATA=as.list(1:25) # Initial DATAs, ?
OM=as.list(1:25) # Initial OMs, ?
XX=NULL
Goal=NULL
NW=5 # Number of weights for model averaging
WM=0.96 # Minimum forgetting weight to use
W=seq(WM,1,length.out=NW) # Weights to use for model averaging
W=round(W,3) # Rounding weights to 3 decimal places
w3=0.999 # Within season Home ground advantage forgetting factor
w_bs=0.8 # Between season forgetting factor
w3_bs=0.999 # Between season Home ground advantage forgetting factor

my.panel.bands <- function(x, y, upper, lower, fill, col,log,
                           subscripts, ..., font, fontface)
{
  upper <- upper[subscripts]
  lower <- lower[subscripts]
  panel.polygon(c(x, rev(x)), c(upper, rev(lower)),
                col = fill, border = FALSE,
                ...)
}

MA_AS=function(ns){
  obj_main=list() # What the function will return
  obj_main$BS=0
  LML=rep(0,NW) # Initial log marginal likelihood for each weight
  Prob=rep(1/NW,NW) # Initial probabilities for each weight value
  season=0
  while(season<ns)
  {
    obj=list() # Object for individual season
    season=season+1
    # Using model averaging for a single
    R=read.csv(NS[season],header=TRUE) # Data file to use
    obj$data=R[,1:10] # Data to use
    Games=dim(obj$data)[1] # Number of games in season
    # Changing $FTR to X, Y and Z, so alphabetically it's Home Win, Draw, Away Win
    for (i in 1:Games){
      if (obj$data$FTR[i]=='H') {
        obj$data$FTR[i] = 'X'
      } else if (obj$data$FTR[i]=='D') {
        obj$data$FTR[i] = 'Y'
      } else {
        obj$data$FTR[i] = 'Z'
      }
    }
    # Renaming of columns 8 and 9. Data to be overwritten with expected H & A goals
    names(obj$data)[8:9]=c("EFTHG","EFTAG")
    obj$data$HomeTeam=factor(obj$data$HomeTeam) # Setting home team column to factors
    obj$data$AwayTeam=factor(obj$data$AwayTeam) # Setting away team column to factors
    obj$data$FTR=factor(obj$data$FTR) # Setting FTR column to factors
    NAMES2=levels(obj$data$HomeTeam) # Names of all teams
    if(season==1){
      NAMES=NAMES2
    }
    obj$data$HomeTeam=as.numeric(obj$data$HomeTeam) # Setting home team column to numeric
    obj$data$AwayTeam=as.numeric(obj$data$AwayTeam) # Setting away team column to numeric
    trunc=4
    ODO[season]=sum(obj$data$FTHG>trunc)+sum(obj$data$FTAG>trunc) # ?
    obj$data$FTHG[obj$data$FTHG>trunc]=trunc # Truncation of home goals to 4 max
    obj$data$FTAG[obj$data$FTAG>trunc]=trunc # Truncation of away goals to 4 max
    obj$data$FTR=as.numeric(obj$data$FTR) # Setting FTR column to numeric: H[1],D[2],A[3]
    obj$data=as.data.frame(obj$data) # Creating a dataframe
    
    A=NULL # Attacking strengths prediction
    AU=NULL # Upper interval of attacking strengths
    AL=NULL # Lower interval of attacking strengths
    B=NULL # Defensive strengths prediction
    BU=NULL # Upper interval of defensive strengths
    BL=NULL # Lower interval of defensive strengths
    G=NULL # Home ground advantage prediction
    GU=NULL # Upper interval of home ground advantage
    GL=NULL # Lower interval of home ground advantage
    PR_HT=rep(0,380) # Initial vector for Pearson residuals of home teams
    PR_AT=rep(0,380) # Initial vector for Pearson residuals of away teams
    PR1=NULL
    PR2=NULL
    
    M=match(NAMES,NAMES2)
    M2=match(NAMES2,NAMES)
    w_surv=which(!is.na(NAMES[M]))
    w_out=which(is.na(NAMES[M]))
    #obj$NAMES[w_surv]
    #obj$NAMES[w_out]
    w_current=which(!is.na(NAMES2[M2]))
    w_in=which(is.na(NAMES2[M2]))
    #obj$NAMES2[w_current]
    #obj$NAMES2[w_in]
    NAMES=NAMES2 # New NAMES is equal to the new list of Teams
    
    if (season>1){
      p_alpha[w_current,]=p_alpha[w_surv,]*w_bs
      p_alpha[w_in,]=rep(32,3)
      q_alpha[w_current,]=q_alpha[w_surv,]*w_bs
      q_alpha[w_in,]=rep(39,3)
      alpha_hat=p_alpha/q_alpha
      
      p_beta[w_current,]=p_beta[w_surv,]*w_bs
      p_beta[w_in,]=rep(39,3)
      q_beta[w_current,]=q_beta[w_surv,]*w_bs
      q_beta[w_in,]=rep(32,3)
      beta_hat=p_beta/q_beta
      
      p_gamma=p_gamma*w3_bs
      q_gamma=q_gamma*w3_bs
      gamma_hat=p_gamma/q_gamma
    } else{
      p_alpha=matrix(10,20,NW) # Initial attacking ps for each weight
      q_alpha=matrix(10,20,NW) # Initial attacking qs for each weight
      alpha_hat=matrix(1,20,NW) # Initial expected attacking strengths for each weight
      p_beta=matrix(10,20,NW) # Initial defensive ps for each weight
      q_beta=matrix(10,20,NW) # Initial defensive qs for each weight
      beta_hat=matrix(1,20,NW) # Initial expected defensive strengths for each weight
      p_gamma=rep(20,NW) # Initial home ground advantage p for each weight
      q_gamma=rep(20,NW) # Initial home ground advantage q for each weight
      gamma_hat=rep(1,NW) # Initial expected home ground advantage for each weight
      LML=rep(0,NW) # Log marginal likelihood for each weight
    }
    
    LML=LML/40
    obj$Omega=NULL
    #Games=max(which(!is.na(obj$data[,4])))
    g=0
    while(g<Games){
      g=g+1 # Next g value for while loop
      D=obj$data[g,] # Game data to use
      HT=D$HomeTeam;AT=D$AwayTeam;X=D$FTHG;Y=D$FTAG
      # Expected score Home
      obj$data[g,8]=rbind(alpha_hat[HT,]*beta_hat[AT,]*gamma_hat) %*%cbind(Prob)
      # Expected score Away
      obj$data[g,9]=rbind(alpha_hat[AT,]*beta_hat[HT,]) %*%cbind(Prob)
      
      #Pearson residuals
      #PR_HT[g]=(X-LHT[g])/sqrt(LHT[g]) # Home Team PR
      #PR_AT[g]=(Y-LAT[g])/sqrt(LAT[g]) # Away Team PR
      
      p_alphaHT= W*p_alpha[HT,] # p_alpha extends for HT
      p_alpha[HT,]=p_alphaHT+rep(X,NW) # p_alpha updates for HT
      q_alphaHT =W*q_alpha[HT,] # q_alpha extends for HT
      # NOTE expectation of beta_hats and gamma_hats are used below
      q_alpha[HT,]=q_alphaHT +beta_hat[AT,]*gamma_hat # q_alpha updates for HT
      alpha_hat[HT,]=p_alpha[HT,]/q_alpha[HT,] # alpha_hat updates for HT
      
      p_alphaAT=W*p_alpha[AT,] # p_alpha extends for AT
      p_alpha[AT,]= p_alphaAT+rep(Y,NW) # p_alpha updates for AT
      q_alphaAT=W*q_alpha[AT,] # q_alpha extends for AT
      # NOTE expectation of beta_hats are used below
      q_alpha[AT,]=q_alphaAT+beta_hat[HT,] # q_alpha updates for HT
      alpha_hat[AT,]=p_alpha[AT,]/q_alpha[AT,] # alpha_hat updates for HT
      
      p_betaHT=W*p_beta[HT,] # p_beta extends for HT
      p_beta[HT,]=p_betaHT+rep(Y,NW) # p_beta updates for HT
      q_betaHT=W*q_beta[HT,] # q_beta extends for HT
      # NOTE expectation of alpha_hats are used below
      q_beta[HT,]=q_betaHT+alpha_hat[AT,] # q_beta updates for HT
      beta_hat[HT,]=p_beta[HT,]/q_beta[HT,] # beta_hat updates for HT
      
      p_betaAT=W*p_beta[AT,] # p_beta extends for AT
      p_beta[AT,]=p_betaAT+rep(X,NW) # p_beta updates for AT
      q_betaAT=W*q_beta[AT,] # q_beta extends for AT
      # NOTE expectation of alpha_hats and gamma_hats are used below
      q_beta[AT,]=q_betaAT+alpha_hat[HT,]*gamma_hat # q_beta updates for AT
      beta_hat[AT,]=p_beta[AT,]/q_beta[AT,] # beta_hat updates for AT
      
      pt_gamma=w3*p_gamma # p_gamma extends
      p_gamma=pt_gamma+rep(X,NW) # p_gamma updates
      qt_gamma=w3*q_gamma # q_gamma extends
      q_gamma=qt_gamma+alpha_hat[HT,]*beta_hat[AT,] # q_gamma updates
      gamma_hat=p_gamma/q_gamma # gamma_hat updates
      
      # Log of marginal mikelihoods updating.
      ML1=lgamma(p_alpha[HT,])-lgamma(p_alphaHT)+p_alphaHT*log(q_alphaHT) +
        -p_alpha[HT,]*log(q_alpha[HT,])
      ML2= lgamma(p_alpha[AT,])-lgamma(p_alphaAT)+p_alphaAT*log(q_alphaAT) +
        -p_alpha[AT,]*log(q_alpha[AT,])
      ML3=lgamma(p_beta[HT,])-lgamma(p_betaHT)+p_betaHT*log(q_betaHT) +
        -p_beta[HT,]*log(q_beta[HT,])
      ML4= lgamma(p_beta[AT,])-lgamma(p_betaAT)+p_betaAT*log(q_betaAT) +
        -p_beta[AT,]*log(q_beta[AT,])
      ML5=lgamma(p_gamma)-lgamma(pt_gamma)+pt_gamma*log(qt_gamma) +
        -p_gamma*log(q_gamma)
      LML=LML+ML1+ML2+ML3+ML4+ML5-lfactorial(X)-lfactorial(Y)
      # representing LMLs as gap to maximum instead (equiv for probabilities)
      LML2=LML-max(LML)
      Prob=exp(LML2)/sum(exp(LML2)) # Probabilities for each W value
      
      if(any(is.nan(Prob))){browser()} # ?
      
      obj$Omega=rbind(obj$Omega,Prob) # Storing of probabilities for each w value
      
      P_alpha=p_alpha %*% cbind(Prob) # p_alpha value to use from MA
      Q_alpha=q_alpha %*% cbind(Prob) # q_alpha value to use from MA
      
      P_beta=p_beta %*% cbind(Prob) # p_beta value to use from MA
      Q_beta=q_beta %*% cbind(Prob) # q_beta value to use from MA
      
      P_gamma=p_gamma %*% cbind(Prob) # p_gamma value to use from MA
      Q_gamma=q_gamma %*% cbind(Prob) # q_gamma value to use from MA
      
      # Confidence intervals appending
      A=rbind(A,qgamma(.5,P_alpha,Q_alpha))
      AU=rbind(AU,qgamma(.75,P_alpha,Q_alpha))
      AL=rbind(AL,qgamma(.25,P_alpha,Q_alpha))
      B=rbind(B,1/qgamma(.5,P_beta,Q_beta))
      BU=rbind(BU,1/qgamma(.75,P_beta,Q_beta))
      BL=rbind(BL,1/qgamma(.25,P_beta,Q_beta))
      G=c(G,qgamma(.5,P_gamma,Q_gamma))
      GU=c(GU,qgamma(.75,P_gamma,Q_gamma))
      GL=c(GL,qgamma(.25,P_gamma,Q_gamma))
    }
    
    Z=I[obj$data$FTR,]
    P=matrix(0,Games,3) # Probability matrix
    P[,2]=dskellam(rep(0,Games),obj$data$EFTHG,obj$data$EFTAG) # prob of draws
    P[,3]=pskellam(rep(-1,Games),obj$data$EFTHG,obj$data$EFTAG) # prob of away wins
    P[,1]=1-P[,2]-P[,3] # prob of home wins
    obj$data$BS=apply((Z-P)^2,1,sum) # Brier score
    obj_main$BS=obj_main$BS+mean(obj$data$BS)
    obj$LS=-apply(Z*log(P) +(1-Z)*log(1-P),1,sum) # Log score
    obj$RPS=RPS(Z,P) # Ranked probability score
    
    obj$Omega=as.data.frame(obj$Omega) # Data frame creation for Omega
    names(obj$Omega)=as.character(W) # Naming each column as the w values
    obj$Omega$season=rep(SEA[season],Games) # Storing of season name
    obj$Omega$game=1:Games # Storing of games numbers
    obj$date=as.Date(obj$data[,2],"%d/%m/%Y") # Storing date column in correct form
    #stack data - Convert from wide to long, for plotting
    AS=stack(as.data.frame(A))[,1]
    AUS=stack(as.data.frame(AU))[,1]
    ALS=stack(as.data.frame(AL))[,1]
    BS=stack(as.data.frame(B))[,1]
    BUS=stack(as.data.frame(BU))[,1]
    BLS=stack(as.data.frame(BL))[,1]
    Club=as.factor(rep(NAMES,Games))
    #levels(Club)=NAMES
    Game=rep(1:Games,each=20) # Used instead of dates
    Date=rep(obj$date,each=20)
    
    # ?
    obj$ALL=data.frame(Game,Club=Club,AS=AS,AUS=AUS,ALS=ALS,BS=BS,BUS=BUS,BLS=BLS,
                       SEA=as.factor(SEA[season])) 
    Game2=1:Games
    obj$ALLG=data.frame(Date,Game=Game2,G=G,GU=GU,GL=GL,
                        SEA=as.factor(SEA[season]))
    
    obj_main$RPS[[season]]=mean(obj$RPS)
    obj_main$LS[[season]]=mean(obj$LS)
    obj_main$All[[season]]=obj$ALL
    obj_main$AllG[[season]]=obj$ALLG
    obj_main$OM[[season]]=obj$Omega
    obj_main$DATA[[season]]=obj$data
    obj_main$NAMESIN[[season]]=NAMES
  }
  return(obj_main)
}

# 2006 to 2020 with model averaging
AS_MA=MA_AS(15)

################################################################################
# Models comparisons - remember to change for each

NW=1 # Number of weights for model averaging
WM=1 # Minimum forgetting weight to use
W=seq(WM,1,length.out=NW) # Weights to use for model averaging
W=round(W,3) # Rounding weights to 3 decimal places
w3=1 # Within season Home ground advantage forgetting factor
w_bs=1 # Between season forgetting factor
w3_bs=1 # Between season Home ground advantage forgetting factor
AS_MA_Static=MA_AS(15)

NW=1 # Number of weights for model averaging
WM=0.9986 # Minimum forgetting weight to use
W=seq(WM,1,length.out=NW) # Weights to use for model averaging
W=round(W,3) # Rounding weights to 3 decimal places
w3=0.9983 # Within season Home ground advantage forgetting factor
w_bs=0.5740 # Between season forgetting factor
w3_bs=0.9983 # Between season Home ground advantage forgetting factor
AS_MA_One=MA_AS(15)

NW=5 # Number of weights for model averaging
WM=0.96 # Minimum forgetting weight to use
W=seq(WM,1,length.out=NW) # Weights to use for model averaging
W=round(W,3) # Rounding weights to 3 decimal places
w3=0.999 # Within season Home ground advantage forgetting factor
w_bs=0.8 # Between season forgetting factor
w3_bs=0.999 # Between season Home ground advantage forgetting factor
AS_MA=MA_AS(15)

ALL_LS=rep(1,15*3)
for (i in 1:15){
  ALL_LS[i]=AS_MA_Static$LS[[i]][1]
}

for (i in 16:30){
  ALL_LS[i]=AS_MA_One$LS[[i-15]][1]
}
for (i in 31:45){
  ALL_LS[i]=AS_MA$LS[[i-15*2]][1]
}

NAMES=c(rep("Static", 15), rep("Single W", 15), rep("Model Averaging", 15))
Dat=data.frame(Mod=NAMES, Mods_LS=ALL_LS, Season=c(SEA, SEA, SEA))
# Dat[Dat$Season == 2006, ]

barchart(Mods_LS~Mod|Season,data=Dat,ylab="Log Scores", scales=list(y=list(rot=0), x=list(rot=90)),col=1:9,main="Log Scores Of All Bayesian Models",as.table=TRUE,grrid=TRUE)

ALL_RPS=rep(1,15*3)
for (i in 1:15){
  ALL_RPS[i]=AS_MA_Static$RPS[[i]][1]
}

for (i in 16:30){
  ALL_RPS[i]=AS_MA_One$RPS[[i-15]][1]
}
for (i in 31:45){
  ALL_RPS[i]=AS_MA$RPS[[i-15*2]][1]
}

NAMES=c(rep("Static", 15), rep("Single Omega_w", 15), rep("Model Averaging", 15))
Dat2=data.frame(Mod=NAMES, Mods_RPS=ALL_RPS, Season=c(SEA, SEA, SEA))
# Dat[Dat$Season == 2006, ]

barchart(Mods_RPS~Mod|Season,data=Dat2, xlab=list("Model", cex=2), ylab=list("Ranked Probability Score", cex=2),
         scales=list(y=list(rot=0), x=list(rot=90)),
         col=1:9,main=list("Ranked Probability Score Of All Bayesian Models", cex=2),as.table=TRUE,grrid=TRUE)

################################################################################

OMEGA=NULL
for(season in 1:15){
  OMEGA=rbind(OMEGA,AS_MA$OM[[season]])
}

xyplot(OMEGA$"1"+OMEGA$"0.99"+OMEGA$"0.98"+OMEGA$"0.97"+OMEGA$"0.96"
       ~game|season,col=1:5,lty=1:5,grid=TRUE,xlab=list("Game Number In Each Season", cex=2), ylab=list("Weights", cex=2),
       main=list("Posterior Weights For The Within Season Forgetting Factor", cex=2),layout=c(8,2),
       data=OMEGA,type="l",auto.key=TRUE,lwd=2,as.table=TRUE,ylim=c(0,1.01),
       key = list(title="Within season forgetting factor",text = list(as.character(W)), 
                  lines= list(lwd=rep(2,NW),lty=NW:1,col=NW:1),columns=3))

LLL=NULL
for(I in 1:15){
  LL=levels(AS_MA$All[[I]][,2])
  LLL=c(LLL,LL)
}
LLL=factor(LLL)
ANAMES=levels(LLL)

################################################################

#Extract All Liverpool scores

#This code plots  the history of AS+DS of each club.

PIT=function(C)
{
 Club=ANAMES[C]
All2=NULL

for(i in 1:15){
  W=which(AS_MA$All[[i]][,2]==Club)
  All2=rbind(All2,AS_MA$All[[i]][W,])
}

names(All2)[3]="Attack"
names(All2)[6]="Defense"

SEA2=paste(as.character(as.numeric(SEA)-1),SEA,sep="-")
lattice.options(default.args = list(as.table = TRUE)) 



xyplot(Defense+Attack~Game|SEA,data=All2,
       main=list("Manchester United Strengths From Single Omega_w Model", cex=2),layout=c(8,2),
       # main=list("Liverpool Strengths From Model Averaging Model", cex=2),layout=c(8,2),
       key = list(title="",text = list(c("Defensive strength","Atacking strength")), 
                  lines= list(lwd=rep(2,2),lty=c(2,1),col=c(4,6)),columns=2) ,
       xlab=list("Game Number In Each Season", cex=2), ylab=list("Strength", cex=2),
       upper = c(All2$BUS,All2$AUS), lower = c(All2$BLS,All2$ALS),
       panel = function(x, y, ...){
         panel.superpose(x, y, panel.groups = my.panel.bands, type='l', col='gray',...)
         panel.xyplot(x, y, type='l',lwd=2,cex=0.6, lty=c(2,1),...)
         panel.grid(v = 5, h = 5)
       }  )
}

Side=match("Man United", ANAMES) # Club number to use
# Side=match("Liverpool", ANAMES) # Club number to use
PIT(Side)

################################################################################

# Strengths Of All Teams For One Season

All3 = AS_MA$All[[2]]
names(All3)[3]="Attack"
names(All3)[6]="Defense"

xyplot(Defense+Attack~Game|Club,data=All3,
       main=list("All Teams Strengths From Model Averaging Model 2006/07", cex=2),layout=c(10,2),
       key = list(title="",text = list(c("Defensive strength","Atacking strength")), 
                  lines= list(lwd=rep(2,2),lty=c(2,1),col=c(4,6)),columns=2) ,
       xlab=list("Game Number In Each Season", cex=2), ylab=list("Strength", cex=2),
       upper = c(All3$BUS,All3$AUS), lower = c(All3$BLS,All3$ALS),
       panel = function(x, y, ...){
         panel.superpose(x, y, panel.groups = my.panel.bands, type='l', col='gray',...)
         panel.xyplot(x, y, type='l',lwd=2,cex=0.6, lty=c(2,1),...)
         panel.grid(v = 5, h = 5)
       }  )

################################################################################

xyplot(G~Adate,data=Goal, xlab="Date",ylab="HGA",groups=SEA,
       upper = Goal$GL, lower =Goal$GU ,ylim=c(1.3,1.41),grid=TRUE,
       main ="Home goal advantage in E.P.L. over more than two decades",
       panel = function(x, y, ...){
         panel.superpose(x, y, panel.groups = my.panel.bands, type='l', col='gray',...)
         panel.xyplot(x, y, type='l',lwd=2,cex=0.6, lty=c(2,1),...)
         panel.grid(v = 5, h = 5)
       }  )

MU=c(1993, 1994, 1996, 1997, 1999, 2000, 2001, 2003, 2007, 2008, 2009, 2011, 2013)
CW=c(2005, 2006, 2010, 2015, 2017)
AW=c(1998, 2002, 2004)
MCW=c(2012, 2014, 2018, 2019)
LW=2016
LIV=2020

ALLPRES=NULL;ALLDRES=NULL
corrs=rep(0,25)
corrs2=rep(0,25)
corCI=matrix(0,25,2)
for(season in 1:25){
  ALLPRES=rbind(ALLPRES,PRES[[season]])
  corrs[season]=cor(ALLPRES[,1],ALLPRES[,2],method="spearman")
  corCI[season,]=as.numeric(cor.test(ALLPRES[,1],ALLPRES[,2])[[9]])
  ALLDRES=rbind(ALLDRES,DRES[[season]])
  corrs2[season]=cor(ALLDRES[,1],ALLDRES[,2],method="spearman")
}
x11()
par(las=2)
par(mfrow=c(1,1))
plot(as.numeric(SEA),corrs,ylim=c(0,.15),main="Correlation of Pearson Residuals",ylab="Spearman c.c.",xlab="season")
grid(10,10,lty=3,col=3)
par(mfrow=c(1,1))

Adate=as.Date(DATA[[2]]$Date,"%d/%m/%Y")
for (S in 3:25)
  Adate=c(Adate,as.Date(DATA[[S]]$Date,"%d/%m/%Y"))

Goal[,1]=Adate
par(las=2)
barplot(ODO,names.arg=SEA,main="Games where one side exceeds 5 goals",ylab="Games with excess of 5 goals")

GRPS=function(season){
thing2=DATA2[[season]]$RPS
thing=DATA[[season]]$RPS
  return(cbind(thing,thing2))}

GBS=function(season){
  thing2=DATA2[[season]]$BS
  thing=DATA[[season]]$BS
  return(cbind(thing,thing2))}

GLS=function(season){
  thing2=DATA2[[season]]$BS
  thing=DATA[[season]]$BS
  return(cbind(thing,thing2))}


BR=as.numeric(SEA)
I=0
while(I<25){
I=I+1;G=GRPS(I);BR[I]=mean(G[,2]>G[,1])
}
DF=data.frame(Season=SEA,BS=BR)
x11()
par(las=2)
barplot(DF[,2],names=SEA)
abline(c(.5,0),col=2,lwd=2)
