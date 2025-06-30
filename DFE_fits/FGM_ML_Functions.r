library(distr);# for fast convolution Loglik

#---------- define distributions using dist package -------------
# density of rescaled y = s/|E(s)|
Dy<-function(n,y0){y0-1/n*Chisq(n,n*y0)}

# density of y = s/|E(s)| isotropic FGM + normal error(0,sye)
Dyobs<-function(n,y0,sye){y0-1/n*Chisq(df=n,ncp=n*y0) + Norm(mean=0,sd=sye)}
# density of s isotropic model + normal error(0,se)
Dsobs<-function(n,y0,se,Es){y0*Es-Es/n*Chisq(df=n,ncp=n*y0) + Norm(mean=0,sd=se)}


# ---------- loglik(obs,params) functions  ----------------------
#obs must be a data.frame with
# unscaled s and measurement error: obs$s and obs$SEM 
# rescaled s and measurement error: obs$y =obs$s/|E(s)| and obs$SEMy = obs$SEM/|E(s)|

# rescaled homoscedastic model: approx Ve = mean(Ve_mutation)
LLobs<-function(obs,n,y0){sum(log(d(Dyobs(n,y0,sqrt(mean(obs$SEMy^2))))(obs$y)))}
# rescaled heteroscedastic model
LLobsHet<-function(obs,n,y0){sum((Vectorize(function(y,sem){log(d(Dyobs(n,y0,sem))(y))})(obs$y,obs$SEMy)))}

# unscaled homoscedastic model:
LLsobs<-function(obs,n,y0,Es){sum(log(d(Dsobs(n,y0,sqrt(mean(obs$SEM^2)),Es))(obs$s)))}

# normal limit (homoscedastic model): approx Ve = mean(Ve_mutation)
LLNobs<-function(obs,Vy){sum(log(d(Norm(mean=0,sd=sqrt(Vy+mean(obs$SEMy^2))))(obs$y)))}

# ------------- brute force maximization on a grid ----------------
MLEygrid<-function(obs,ns,y0s){
MLEs=as.data.frame(matrix(0,ncol=3,nrow=length(y0s)*length(ns)))
names(MLEs)=c('n','y0','LogLik')
i=0;
for(n in ns){
  print(n);
  for(y0 in y0s){
    i=i+1;
    MLEs[i,]=c(n,y0,LLobs(obs,n,y0));
  }
}
return(MLEs)
}
