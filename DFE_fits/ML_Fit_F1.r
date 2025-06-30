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
F10 <- read_delim("data_txt_files/F1_single_mutations.txt", 
                                          delim = "\t", escape_double = FALSE, 
                                          trim_ws = TRUE)

nmut=dim(F10)[1] # number of mutations in dataset

# P(lethal) = f(type) (type: synonymous (S), non-synonymous (NS), stop codon (STOP), intergenic (I)
viabilities=ftable(F10[,c('type','viable')])
# fisher exact test of independence of viability vs type 
fisher.test(viabilities) # p <<1 by fisher exact test (but again small classes is still a pb for power)

#do the same for just two types : synonymous (S+I) or not (NS + STOP)
F10$syn=ifelse(F10$type %in% c('S','I'),1,0)
viabilities2=ftable(F10[,c('syn','viable')])
#fisher exact test of independence: syn or not  vs. viability
fisher.test(viabilities2) # p <<1 only the non syn are lethal
# => estimate P(lethal) over all mutations: 19%
Plethal=sum(viabilities2[,1])/sum(viabilities2) 
# CI for P(lethal) based on binomial sample [13%-29%]
c(qbinom(0.025,size=nmut,prob=Plethal),qbinom(0.975,size=nmut,prob=Plethal))/nmut

# now focus on the viable mutations (VSV0->VSV)
F1=F10[F10$viable==1,]
Es=abs(mean(F1$s))
F1$y=F1$s/Es
F1$SEMy=F1$SEM/Es
View(F1)

#---------- explore the dataset (Ve(y) vs |y| estimate)-------------------------
plot(abs(F1$y),F1$SEMy^2,log='x',xlab='|y|',ylab='Ve(y)')
abline(h=mean(F1$SEMy),lty='dashed')
# there is no significant heteroscedasticity found (SEM² vs |s|, loglog)
hetTest=lm(log(F1$SEMy^2) ~ log(abs(F1$y)))
summary(hetTest)
# so no need for heteroscedastic MLE here

#remove 1 OUTLIER with extreme error variance
F1c=F1[-which.max(F1$SEMy),]
lines(abs(F1c$y),F1c$SEMy^2,col='blue',type='p',pch=16)
abline(h=mean(F1c$SEMy^2),lty='dashed',col='blue')
lines(abs(F1[which.max(F1$SEMy),]$y),F1[which.max(F1$SEMy),]$SEMy^2,type='p',pch=16,col='red')

# change rescaling after removing OUTLIER
Esc=abs(mean(F1c$s))
F1c$y=F1c$s/Esc
F1c$SEMy=F1c$SEM/Esc

#----------- MLE exploration for F1 (sanjuan 2004) ----------------------------------------
y0s=seq(from=0,to=10*max(F1c$s),length.out=100)
ns=1:10;
MLEs=MLEygrid(F1c,ns,y0s)  #compute loglik profile

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

#LRT for y0 = 0: significant y0 > 0
X=2*(LLobs(F1c,best$n,best$y0)-LLobs(F1c,best$n,0))
1-pchisq(X,df=1)

# -----------   MLE using maxLik package ---------------------
## estimate mean and variance of normal random vector
set.seed(123)
llf<-function(par){
  logy0<-par[1]
  a<-par[2]
  return(LLsobs(F1c,best$n,0.001+exp(logy0),Esc*exp(a)))
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
write.table(estimates,file="F1_estimates.txt")


#correlation between estimates: PB
cov2cor(vcov(ml))







