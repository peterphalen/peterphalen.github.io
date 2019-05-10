graphics.off()
rm(list=ls(all=TRUE))

library(rstan)

#------------------------------------------------------------------------------
# THE DATA.

# Temperature data from
# http://academic.udayton.edu/kissock/http/Weather/default.htm
dataMat = read.table( "https://www.peterphalen.com/temp/WIMADISO.txt" , 
                      col.names=c("Month","Date","Year","AveTemp") )
cityName = "Madison, Wisconsin"

# Re-code missing data from -99 to NA:
dataMat[ dataMat[,"AveTemp"]==(-99) , "AveTemp" ] = NA

# Plot data:
plot( dataMat$AveTemp , main=cityName ,
      xlab="Day since 1/1/95" , ylab="Ave Daily Temp (F)" , type="l" ,
      xaxt="n")
Jan1RowIdx = which( dataMat$Month==1 & dataMat$Date==1 )
Jul1RowIdx = which( dataMat$Month==7 & dataMat$Date==1 )
for ( i in 1:length(Jan1RowIdx) ) {
  abline( v=Jan1RowIdx[i]-0.5 , col="grey" )
}
for ( i in 1:length(Jul1RowIdx) ) {
  text( Jul1RowIdx[i] , min(dataMat$AveTemp,na.rm=TRUE) , 
        dataMat[Jul1RowIdx[i],"Year"] , adj=c(0.5,0), cex=.5 )
}

# Re-name data for use in Stan model:
x = 1:nrow(dataMat)
y = as.vector(dataMat[,"AveTemp"])

# Clip to integer number of cycles, to minimize influence of end effects:
lastIdx = max( which( dataMat$Month==12 & dataMat$Date==31 ) )
x = x[1:lastIdx]
y = y[1:lastIdx]

# Remove missing points:
includeIdx = which( !is.na( y ) )
x = x[includeIdx]
y = y[includeIdx]

# Because initial day is arbitrary, re-center at mean of x values:
x = round( x - mean(x) )
centerDateIdx = which( x == 0 )
centerMonth = dataMat[centerDateIdx,"Month"]
centerDate = dataMat[centerDateIdx,"Date"]
centerYear = dataMat[centerDateIdx,"Year"]

# # For debugging/testing, include only first few years:
# limit = 4*365
# x = x[1:limit]
# y = y[1:limit]

# Specify data, as a list.
daysPerYear = 365.24219 # Tropical Year. Will be used instead of estimating wl.
dataList = list(
  x = x ,
  y = y ,
  N = length(x),
  wl = daysPerYear/(2*pi) 
)

#------------------------------------------------------------------------------
# THE MODEL.
modelstring = "
data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
  real wl;
}

parameters{
  real<lower=-1.2,upper=1.1> ar1;
  real<lower=0> sigma;
  real<lower=0> amp;
  real<lower=-183,upper=183> phase; // plus and minus half cycle
  real beta0;
  real beta1;
  real<lower=0> nu;
//  real wl;  // commented out
}

model {
  vector[N] trend; 
  vector[N] mu;

  trend[1] = beta0 + beta1 * x[1] + amp * cos( ( x[1] - phase ) / wl );
  mu[1] = trend[1];
  
  trend[2:N] = beta0 + beta1 * x[2:N] + amp * cos( ( x[2:N] - phase ) / wl );
  mu[2:N] = trend[2:N] + ar1 * ( y[1:(N-1)] - trend[1:(N-1)] );
    
  ar1 ~ normal( 0 , 0.5 );
  beta0 ~ normal( 45 , 10 );
  beta1 ~ normal( 0 , 2 );
  sigma ~ normal( 0 , 20 );
  amp ~ normal( 25 , 10 );
  phase ~ normal( 0 , 50 ); 
  nu ~ exponential(0.04);
 // wl ~ normal(58.1,.25);
  y ~ student_t(nu, // df
                mu, 
                sigma);
}
" # close quote for modelstring

fit <- stan(model_code=modelstring, data=dataList, 
            cores=4, iter=2000)

#------------------------------------------------------------------------------
# EXAMINE THE RESULTS

mcmcChain = as.data.frame( fit )

elapsed_time <- round(max(rowSums(get_elapsed_time(fit)))/60,1)
elapsed_time <- paste("model took",elapsed_time,
                      "minutes to fit in stan")

# Plot data with posterior predictive curves:
layout( matrix( c(rep(1,4),1+1:8) , nrow=3 , byrow=TRUE) , heights=c(2,1) )
plot( x,y , main=cityName , cex.lab=1.4 ,
      xlab=elapsed_time , 
      ylab="Ave. Daily Temp. (deg. F)" , type="l" , xaxt="n")
Jan1RowIdx = which( dataMat$Month==1 & dataMat$Date==1 )
Jul1RowIdx = which( dataMat$Month==7 & dataMat$Date==1 )
for ( i in 1:length(Jan1RowIdx) ) {
  abline( v=Jan1RowIdx[i]+x[1]-0.5 , col="grey" )
}
for ( i in 1:length(Jul1RowIdx) ) {
  text( Jul1RowIdx[i]+x[1] , min(dataMat$AveTemp,na.rm=TRUE) , 
        dataMat[Jul1RowIdx[i],"Year"] , adj=c(0.5,0) )
}
curvesToPlot = round(seq(1,nrow(mcmcChain),length=20))
for ( i in curvesToPlot ) {
 lines( x , 
        mcmcChain[i,"beta0"] + mcmcChain[i,"beta1"]*x + 
          mcmcChain[i,"amp"] * cos( ( x - mcmcChain[i,"phase"] ) / dataList$wl ) ,
        col="skyblue" , lwd=2 )
}
hist( mcmcChain[,"beta0"] , xlab="Deg. F" , main="Intercept")       
hist( mcmcChain[,"beta1"]*daysPerYear , xlab="Deg. F / Year" ,
                     main="Linear Trend"  )
hist( mcmcChain[,"amp"] , xlab="Deg. F" , main="Amplitude" )
hist( mcmcChain[,"phase"] , 
                     xlab=paste("Days since",centerMonth,"/",centerDate) , 
                     main="Peak Temp. Day" )
hist( mcmcChain[,"ar1"] , xlab="AR(1) Coef." , main="AR(1) Coef." )
hist( mcmcChain[,"sigma"] , xlab="Deg. F" , main="SD noise" )
hist( mcmcChain[,"nu"] , xlab="student t df" , main="Normality"  )
# hist( mcmcChain[,"wl"] * 2 * pi, xlab="days" , main="empirical year length" )

#------------------------------------------------------------------------------
