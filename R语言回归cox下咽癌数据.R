library(foreign)
library(rms)
mydata<-read.spss("C:/Users/CSH/Desktop/LAHPSCC.sav")
mydata<-as.data.frame(mydata)
mydata<-na.omit(mydata)
head(mydata)
mydata$Status<-ifelse(mydata$Status=="Dead",1,0)
mydata$CSS<-ifelse(mydata$CSS=="Dead",1,0)
mydata$Race<-relevel(mydata$Race,ref = 'white')
mydata$Race<-relevel(mydata$Sex,ref = 'female')
mydata$Race<-relevel(mydata$Grade,ref = 'I-II')
mydata$Site<-relevel(mydata$Site,ref = 'Pyriform sinus')
mydata$Tstage<-relevel(mydata$Tstage,ref = 'T1')
mydata$N<-relevel(mydata$N,ref = 'N0')
mydata$Treatment<-relevel(mydata$Treatment,ref = 'S')
mydata$clinical<-relevel(mydata$clinical,ref = 'III')
dd<-datadist(mydata)
options(datadist = 'dd')
coxm<-cph(Surv(times,CSS==1)~Age+Race+Site+Sex+Tstage+N+Treatment+clinical,x=T,y=T,data = mydata,surv = T)
coxm1<-cph(Surv(times,CSS==1)~Tstage+N,x=T,y=T,data = mydata,surv = T)
coxm2<-cph(Surv(times,CSS==1)~clinical,x=T,y=T,data = mydata,surv = T)
library(Hmisc)
surv<-Survival(coxm)
rcorrcens(surv(times,CSS)~predict(coxm),data = mydata)
surv<-Survival(coxm1)
rcorrcens(surv(times,CSS)~predict(coxm1),data = mydata)
surv<-Survival(coxm2)
rcorrcens(surv(times,CSS)~predict(coxm2),data = mydata)

surv1<-function(x)surv(1*12,lp=x)#defined time.inc,1 year OS
surv2<-function(x)surv(1*36,lp=x)#defined time.inc,3 year OS
surv3<-function(x)surv(1*60,lp=x)#defined time.inc,5 year OS
nom<-nomogram(coxm,fun = list(surv1,surv2,surv3),lp = F,funlabel = c("1-Year OS",'3-Year OS','5-YearOS'),maxscale = 100,fun.at = c('0.95','0.85','0.80','0.70','0.60','0.50','0.4','0.3','0.2','0.1'))
plot((nom),xfrac=.3)

cal<-calibrate(coxm,cmethod = 'KM',method = 'boot',u=60,m=150,B=1000)
plot(cal,lwd=2,lty=1,errbar.col=c(rgb(0,118,192,maxColorValue = 255)),
     xlim=c(0.0,1),ylim=c(0.0,1),
     xlab="Nomogram-predicted Probability of 5-year OS",
     ylab="Actual 5-year OS(proportion)",
     col=c(rgb(192,98,83,maxColorValue = 255)))
lines(cal[,c("mean.predicted","KM")],type = "b",lwd=2,col=c(rgb(192,98,83,maxColorValue = 255)),pch=16) 
abline(0,1,lty=3,lwd=2,col=c(rgb(0,118,192,maxColorValue = 255)))


library('nomogramFormula')
ret = nomogramFormula::TotalPoints.rms(rd=mydata,fit=coxm,nom=nom)
ret$'total points'
write.csv(ret,'LAHPSCC.csv')



library(timeROC)
mydata$model1 = predict(coxm,data = mydata)
mydata$model2 = predict(coxm1,data = mydata)

roc1 = timeROC(T=mydata$times,
               delta = mydata$Status,
               marker = mydata$model1,
               cause = 1,
               times = 60,
               ROC =TRUE,
               iid = TRUE
)
roc1$AUC
roc2 = timeROC(T=mydata$times,
               delta = mydata$Status,
               marker = mydata$model2,
               cause = 1,
               times = 60,
               ROC =TRUE,
               iid = TRUE
)
roc2$AUC     
pdf('ROC.pdf')
plot(roc1,time=60,col='red',title ='ROC',lw=2)
plot(roc2,time=60,add=T,col='green',lw=2)
legend('bottomright',
       c(paste('model1 AUC',round(roc1$AUC,3)[2],sep=':'),
         paste('model2 AUC',round(roc2$AUC,3)[2],sep=':'))
       ,lw=2,col=c('red','green'))

dev.off()
compare(roc1,roc2)



library(nricens)
library(survival)
library(foreign)
mydata<-read.spss("C:/Users/CSH/Desktop/LAHPSCC.sav")
mydata<-as.data.frame(mydata)
mydata<-na.omit(mydata)
mydata$Status<-ifelse(mydata$Status=="Dead",1,0)
names(mydata)
time<-mydata$Months
status<-mydata$Status
j1<-as.matrix(subset(mydata,select = c(Tstage,N,M)))
j2<-as.matrix(subset(mydata,select = c(Age,Race,Site,Tstage,N,M,Treatment)))
mod.std<-coxph(Surv(time,status)~.,data.frame(time,status,j1),x=TRUE)
mod.new<-coxph(Surv(time,status)~.,data.frame(time,status,j2),x=TRUE)
p.std=get.risk.coxph(mod.std,t0=60)
p.new=get.risk.coxph(mod.new,t0=60)
nricens(mdl.std = mod.std,mdl.new = mod.new,t0=60,cut = c(0.2,0.4),
        niter = 100,alpha = 0.05,updown = 'category')


library(nricens)
library(survival)
library(foreign)
mydata<-read.spss("C:/Users/CSH/Desktop/LAHPSCC.sav")
dat=mydata[1:1117]
dat$time=as.numeric(dat$times)
dat$status=ifelse(dat$Status=="Dead",1,0)
t0=12*5

indata0=as.matrix(subset(dat,select = c(clinical)))
indata1=as.matrix(subset(dat,select = c(Age,Race,Site,Tstage,N,M,Treatment)))
covs0<-as.matrix(indata0[,c(-1,-2)])
covs1<-as.matrix(indata1[,c(-1,-2)])
library(survIDINRI)
x<-IDI.INF(mydata[,2:1],covs0,covs1,t0,npert = 100)


library(dcurves)
library(foreign)
library(survival)
library(dplyr)
mydata<-read.spss("C:/Users/CSH/Desktop/LAHPSCC.sav",use.value.labels = F,to.data.frame = T)
mydata<-na.omit(mydata)
mydata$Age<-as.factor(mydata$Age)
mydata$Race<-as.factor(mydata$Race)
mydata$Site<-as.factor(mydata$Site)
mydata$Tstage<-as.factor(mydata$Tstage)
mydata$N<-as.factor(mydata$N)
mydata$M<-as.factor(mydata$M)
mydata$Treatment<-as.factor(mydata$Treatment)
f1<-coxph(Surv(Months,Status)~Tstage+N+M,mydata)
f2<-coxph(Surv(Months,Status)~Tstage+N+M+Age+Race+Site+Treatment,mydata)
mydata$pr_failuref136=c(1-(summary(survfit(f1,newdata=mydata),Months=36)$surv))
mydata$pf236=c(1-(summary(survfit(f1,newdata=mydata),Months=36)$surv))
mydata$pf336=c(1-(summary(survfit(f1,newdata=mydata),Months=36)$surv))



library(rms)
library(ggDCA)
library(survival)  
library(foreign)
rm(list = ls())     
mydata<-read.spss("C:/Users/CSH/Desktop/LAHPSCC.sav")
mydata<-as.data.frame(mydata)
head(mydata)
attach(mydata)
dd<-datadist(mydata)
options(datadist='dd')
model1<-coxph(Surv(Months,Status==1)~Tstage+N+M+Age+Race+Site+Treatment,data=mydata)        
dca1<-dca(model1,
          new.data = NULL,
          Months=60)      



library(rms)
library(foreign)
library(ggDCA)
library(ggplot2)
mydata<-read.spss("C:/Users/CSH/Desktop/LAHPSCC.sav")
mydata<-as.data.frame(mydata)
mydata<-na.omit(mydata)
head(mydata)
mydata$Status<-ifelse(mydata$Status=="Dead",1,0)
mydata$Race<-relevel(mydata$Race,ref = 'white')
mydata$Site<-relevel(mydata$Site,ref = 'Pyriform sinus')
mydata$Tstage<-relevel(mydata$Tstage,ref = 'T1')
mydata$N<-relevel(mydata$N,ref = 'N0')
mydata$M<-relevel(mydata$M,ref = 'M0')
mydata$Treatment<-relevel(mydata$Treatment,ref = 'S')
mydata$clinical<-relevel(mydata$clinical,ref = 'III')
dd<-datadist(mydata)
options(datadist = 'dd')
M1<-cph(Surv(times,Status)~clinical,data=mydata)
M2<-cph(Surv(times,Status)~Age+Site+Race+Tstage+N+Treatment,data=mydata)
d<-dca(M1,M2,
       times=c(1))

ggplot(d,linetype=1)  

