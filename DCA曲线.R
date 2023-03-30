setwd("C:/Users/CSH/Desktop")
source("stdca.R") #一定要把stdca.R放在之前设定的起始目录中。
dev<-read.csv("data.csv")
head(dev)
str(dev)
library(rms)
library(foreign)
library(survival)

Srv=Surv(dev$time,dev$Status)
coxmod=coxph(Srv~Tstage+N+Clinical+EBV+VGTV+VLN+GTVADCDMFS+LNADCDMFS,data=dev)
dev$three.years.Survival.Probabilitynew = c(1- (summary(survfit(coxmod,newdata=dev),times=24)$surv))  #计算3年生存率
dev$five.years.Survival.Probabilitynew = c(1- (summary(survfit(coxmod,newdata=dev),times=12)$surv))  #计算5年生存率
write.csv(dev, "devnnew.csv") 
stdca(data=dev,outcome="Status",ttoutcome="time",timepoint=12,predictors="three.years.Survival.Probabilitynew",xstop=0.8,smooth=TRUE)
stdca(data=dev,outcome="Status",ttoutcome="time",timepoint=24,predictors="five.years.Survival.Probabilitynew",xstop=0.8,smooth=TRUE)


#两模型比较
coxmod1<-coxph(Srv ~Tstage+N+Clinical+EBV+VGTV+VLN+GTVADCDMFS+LNADCDMFS,data=dev)
coxmod2<-coxph(Srv ~Tstage+N,data=dev)

dev$model1<-c(1- (summary(survfit(coxmod1,newdata=dev),times=24)$surv))

dev$model2<-c(1- (summary(survfit(coxmod2,newdata=dev),times=24)$surv))



stdca(data=dev,outcome="Status",ttoutcome="time",timepoint=24,predictors=c("model1"),smooth=TRUE)
stdca(data=dev,outcome="Status",ttoutcome="time",timepoint=24,predictors=c("model1","model2"),xstop = 0.8,smooth=TRUE)
