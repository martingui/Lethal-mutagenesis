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
VSV0 <- read_delim("data_txt_files/VSV_random_single_mutations.txt", 
                   delim = "\t", escape_double = FALSE, 
                   trim_ws = TRUE)

nmut=dim(VSV0)[1] # number of mutations in dataset

# P(lethal) = f(type) (type: synonymous (S), non-synonymous (NS), stop codon (STOP), intergenic (I)
viabilities=ftable(VSV0[,c('type','viable')])
# fisher exact test of independence of viability vs type 
fisher.test(viabilities) # p <<1 by fisher exact test (but again small classes is still a pb for power)

#do the same for just two types : synonymous (S+I) or not (NS + STOP)
VSV0$syn=ifelse(VSV0$type %in% c('S','I'),1,0)
viabilities2=ftable(VSV0[,c('syn','viable')])
#fisher exact test of independence: syn or not  vs. viability
fisher.test(viabilities2) # p <<1 only the non syn are lethal
# => estimate P(lethal) over all mutations: 19%
Plethal=sum(viabilities2[,1])/sum(viabilities2) 
# CI for P(lethal) based on binomial sample [13%-29%]
c(qbinom(0.025,size=nmut,prob=Plethal),qbinom(0.975,size=nmut,prob=Plethal))/nmut

# now focus on the viable mutations (VSV0->VSV)
VSV=VSV0[VSV0$viable==1,]
Es=abs(mean(VSV$s))
VSV$y=VSV$s/Es
VSV$SEMy=VSV$SEM/Es
View(VSV)

#---------- explore the dataset (Ve(y) vs |y| estimate)-------------------------
plot(abs(VSV$y),VSV$SEMy^2,log='x',xlab='|y|',ylab='Ve(y)')
abline(h=mean(VSV$SEMy),lty='dashed')
# there is no significant heteroscedasticity found (SEM² vs |s|, loglog)
hetTest=lm(log(VSV$SEMy^2) ~ log(abs(VSV$y)))
summary(hetTest)
# so no need for heteroscedastic MLE here

#remove 1 OUTLIER with extreme error variance
VSVc=VSV[-which.max(VSV$SEMy),]
lines(abs(VSVc$y),VSVc$SEMy^2,col='blue',type='p',pch=16)
abline(h=mean(VSVc$SEMy^2),lty='dashed',col='blue')
lines(abs(VSV[which.max(VSV$SEMy),]$y),VSV[which.max(VSV$SEMy),]$SEMy^2,type='p',pch=16,col='red')

# change rescaling after removing OUTLIER
Esc=abs(mean(VSVc$s))
VSVc$y=VSVc$s/Esc
VSVc$SEMy=VSVc$SEM/Esc

#----------- MLE exploration for VSV (sanjuan 2004) ----------------------------------------
y0s=seq(from=0,to=10*max(VSVc$s),length.out=100)
ns=1:10;
MLEs=MLEygrid(VSVc,ns,y0s)  #compute loglik profile

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

ggsave("VSV_LL_profiles.pdf")

#LRT for y0 = 0: marginally significant y0 > 0 p = 0.088
X=2*(LLobs(VSVc,best$n,best$y0)-LLobs(VSVc,best$n,0))
1-pchisq(X,df=1)

# -----------   MLE using maxLik package ---------------------
## estimate mean and variance of normal random vector
set.seed(123)
llf<-function(par){
  logy0<-par[1]
  a<-par[2]
  return(LLsobs(VSVc,best$n,0.001+exp(logy0),Esc*exp(a)))
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
write.table(estimates,file="VSV_estimates.txt")


#correlation between estimates: PB
cov2cor(vcov(ml))







