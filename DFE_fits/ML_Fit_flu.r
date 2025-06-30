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


#----------- import data for flu (visher 2016) --------------------------------
#import flu dataset (caution with . vs , as numeric separator)
flu0 <- read_delim("data_txt_files/flu.csv", 
                   delim = ",", escape_double = FALSE, 
                   trim_ws = TRUE)

colnames(flu0) = c('mut','s','SEM')

flu0$s=flu0$s-1

nmut=dim(flu0)[1] # number of mutations in dataset

# now focus on the viable mutations (flu0->flu)

flu=flu0
Es=abs(mean(flu$s))
flu$y=flu$s/Es
flu$SEMy=flu$SEM/Es
View(flu)

Plethal=0.316

#---------- explore the dataset (Ve(y) vs |y| estimate)-------------------------
plot(abs(flu$y),flu$SEMy^2,log='x',xlab='|y|',ylab='Ve(y)')
abline(h=mean(flu$SEMy),lty='dashed')
# there is no significant heteroscedasticity found (SEM² vs |s|, loglog)
#hetTest=lm(log(flu$SEMy^2) ~ log(abs(flu$y)))
#summary(hetTest)
# so no need for heteroscedastic MLE here


#----------- MLE exploration for flu (visher 2016) ----------------------------------------
y0s=seq(from=0,to=10*max(flu$s),length.out=100)
ns=1:10;
MLEs=MLEygrid(flu,ns,y0s)  #compute loglik profile

#find the lowest -loglik
best=MLEs[which.max(MLEs$LogLik),]
best


# plot LL profiles across n integer values and y0
best=MLEs[which.max(MLEs$LogLik),]
par(mfrow=c(1,1))

ggplot(MLEs, aes(x=y0 ,y = LogLik, color = factor(n))) +
  geom_line() +
  theme_minimal() +
  annotate("point", x = best$y0, y = best$LogLik, size = 3, color = "red") +
  labs(title = "logLik vs maladaptation",
       x = "y0",
       y = "logLik")

ggsave("flu_LL_profiles.pdf")

#LRT for y0 = 0: marginally significant y0 > 0 p = 0.088
X=2*(LLobs(flu,best$n,best$y0)-LLobs(flu,best$n,0))
1-pchisq(X,df=1)

# -----------   MLE using maxLik package ---------------------
## estimate mean and variance of normal random vector
set.seed(123)
llf<-function(par){
  logy0<-par[1]
  a<-par[2]
  return(LLsobs(flu,best$n,0.001+exp(logy0),Esc*exp(a)))
  }
# check LL definition at guess
llf(c(log(best$y0),0))
# perform MLE
ml <- maxLik(llf, start = c(logy0=log(best$y0),a=0),method="NM")
# summary of the ML fit
summary(ml)
# estimated parameter values and their CI
estimates=data.frame(what=c("MLE","2.5%","97.5%"),n=best$n,y0=best$y0,Es=Esc,Plethal=Plethal)
estimates[1,3:4]=c(exp(coef(ml)[1])+0.001,Esc*exp(coef(ml)[2]))
estimates[2:3,3]=exp(confint(ml)[1,])+0.001
estimates[2:3,4]=Esc*exp(confint(ml)[2,])
estimates[2:3,5]=c(qbinom(0.025,size=nmut,prob=Plethal),qbinom(0.975,size=nmut,prob=Plethal))/nmut
View(estimates)
write.table(estimates,file="flu_estimates.txt")


#correlation between estimates: PB
cov2cor(vcov(ml))







