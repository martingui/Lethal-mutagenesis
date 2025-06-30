#load packages
library(stats4);library(MASS);library(graphics);
library(readr) # to read .txt files
#library(rgenoud);# for genetic algo maximization
library(distr);# for fast convolution Loglik
library(ggplot2)
library(gridExtra) # to show multiple ggplots together
library(tictoc)  # to measure CPU time
library(maxLik) # for mle diagnostics

# request source functions
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # get current script directory
#source functions that are in the same directory
source("FGM_ML_Functions.r")


#----------- import data for VSV (sanjuan 2004) --------------------------------
#import VSV dataset (caution with . vs , as numeric separator)
PHIX0 <- read_delim("data_txt_files/PHIX174_single_mutations.txt", 
                                          delim = "\t", escape_double = FALSE, 
                                          trim_ws = TRUE)

nmut=dim(PHIX0)[1] # number of mutations in dataset

# P(lethal) = f(type) (type: synonymous (S), non-synonymous (NS), stop codon (STOP), intergenic (I)
viabilities=ftable(PHIX0[,c('type','viable')])
# dependence of type and viability 
fisher.test(viabilities)

# => estimate P(lethal) over all mutations: 20%
Plethal=sum(viabilities[,1])/sum(viabilities) 
# CI for P(lethal) based on binomial sample [9%-31%]
c(qbinom(0.025,size=nmut,prob=Plethal)/nmut,Plethal,qbinom(0.975,size=nmut,prob=Plethal)/nmut)

# now focus on the viable mutations (VSV0->VSV)
PHIX=PHIX0[PHIX0$viable==1,]
Es=abs(mean(PHIX$s))
PHIX$y=PHIX$s/Es
PHIX$SEMy=PHIX$SEM/Es
View(PHIX)

#---------- explore the dataset (Ve(y) vs |y| estimate)-------------------------
plot(abs(PHIX$y),PHIX$SEMy^2,log='',xlab='|y|',ylab='Ve(y)',pch=16,col='red')
abline(h=mean(PHIX$SEMy^2),lty='dashed')
# no significant heteroscedasticity found (SEM² vs |s|)
hetTest=lm(PHIX$SEMy^2 ~ abs(PHIX$y))
summary(hetTest)
# so no need for heteroscedastic MLE here

#no obvious OUTLIERS with extreme error variance
PHIXc=PHIX#[order(-PHIX$SEMy),][-(1:2),]
lines(abs(PHIXc$y),PHIXc$SEMy^2,col='blue',type='p',pch=16)
abline(h=mean(PHIXc$SEMy^2),lty='dashed',col='blue')

# change rescaling after removing OUTLIER
Esc=abs(mean(PHIXc$s))
PHIXc$y=PHIXc$s/Esc
PHIXc$SEMy=PHIXc$SEM/Esc

hist(PHIXc$y,br=20)

#----------- MLE exploration for PHIX (sanjuan 2004) ----------------------------------------
y0s=seq(from=0,to=20*max(PHIXc$s),length.out=100)
ns=1:10;
MLEs=MLEygrid(PHIXc,ns,y0s)  #compute loglik profile

#find the lowest -loglik
best=MLEs[which.max(MLEs$LogLik),]
best


# plot LL profiles across n integer values and y0
par(mfrow=c(1,1))

ggplot(MLEs, aes(x=y0 ,y = LogLik, color = factor(n))) +
  geom_line() +
  theme_minimal() +
  annotate("point", x = best$y0, y = best$LogLik, size = 3, color = "red") +
  geom_hline(yintercept=best$LogLik-1,
             linetype="dashed", color = "blue")+
  annotate("point", x = best$y0, y = best$LogLik, size = 3, color = "red") +
  labs(title = "logLik vs maladaptation",x = "y0",y = "logLik")

ggsave("PHIX_LL_profiles.pdf")


#LRT for y0 = 0: non significant y0
X=2*(LLobs(PHIXc,best$n,best$y0)-LLobs(PHIXc,best$n,0))
1-pchisq(X,df=1)


# -----------   MLE using maxLik package ---------------------
# maximize for {y0,Es} at n = 1
set.seed(123)
llf<-function(par){
  logy0<-par[1]
  a<-par[2]
  return(LLsobs(PHIXc,best$n,1E-8+exp(logy0),Esc*exp(a)))
  }

# check LL definition at guess
llf(c(log(best$y0),0))
# perform MLE
ml <- maxLik(llf, start = c(logy0=log(best$y0),a=0),method="NM")
# summary of the ML fit: : pb with std errors
summary(ml)


# but we have seen that y0 is not significant a priori 
# maximize for {n,Es} at y0 = 0 for n = 1 or more
set.seed(123)
llf<-function(par){
  logn<-par[1]
  a<-par[2]
  return(LLsobs(PHIXc,1.001+exp(logn),0,Esc*exp(a)))
}

# check LL definition at guess
llf(c(-10,0))
# perform MLE
ml1 <- maxLik(llf, start = c(logn=-10,a=0),method="NM")
# summary of the ML fit: OK proper convergence and std errors
summary(ml1)



# estimated parameter values and their CI
estimates=data.frame(what=c("MLE","2.5%","97.5%"),n=best$n,y0=0,Es=Esc,Plethal=Plethal)
estimates[1,c(2,4)]=c(1.001+exp(coef(ml1)[1]),Esc*exp(coef(ml1)[2]))
estimates[2:3,2]=1.001+exp(confint(ml1)[1,])
estimates[2:3,4]=Esc*exp(confint(ml1)[2,])
estimates[2:3,5]=c(qbinom(0.025,size=nmut,prob=Plethal),qbinom(0.975,size=nmut,prob=Plethal))/nmut
View(estimates)
write.table(estimates,file="PHIX_estimates.txt")






