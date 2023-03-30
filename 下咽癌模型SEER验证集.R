library(foreign)
library(rms)
mydata<-read.spss("C:/Users/CSH/Desktop/LAHPSCC.sav")
mydata<-as.data.frame(mydata)
mydata<-na.omit(mydata)
head(mydata)
mydata$status<-ifelse(mydata$status=="yes",1,0)
mydata$risk<-relevel(mydata$risk,ref = 'low')
mydata$Site<-relevel(mydata$Site,ref = 'Pyriform sinus')
mydata$Tstage<-relevel(mydata$Tstage,ref = 'T1')
mydata$N<-relevel(mydata$N,ref = 'N0')
mydata$M<-relevel(mydata$M,ref = 'M0')
mydata$Treatment<-relevel(mydata$Treatment,ref = 'S+R+C')
dd<-datadist(mydata)
options(datadist = 'dd')
coxm<-cph(Surv(Months,Status==1)~Age+Race+Site+Tstage+N+M+Treatment,x=T,y=T,data = mydata,surv = T)
coxm1<-cph(Surv(Months,Status==1)~Tstage+N+M,x=T,y=T,data = mydata,surv = T)
library(Hmisc)
surv<-Survival(coxm)
surv<-Survival(coxm1)
rcorrcens(surv(Months,Status)~predict(coxm),data = mydata)
rcorrcens(surv(Months,Status)~predict(coxm1),data = mydata)
coxm1
coxm

surv1<-function(x)surv(1*12,lp=x)#defined time.inc,1 year OS
surv2<-function(x)surv(1*36,lp=x)#defined time.inc,3 year OS
surv3<-function(x)surv(1*60,lp=x)#defined time.inc,5 year OS
nom<-nomogram(coxm,fun = list(surv1,surv2,surv3),lp = F,funlabel = c("1-Year OS",'3-Year OS','5-YearOS'),maxscale = 100,fun.at = c('0.95','0.85','0.80','0.70','0.60','0.50','0.4','0.3','0.2','0.1'))
plot((nom),xfrac=.3)

cal<-calibrate(coxm,cmethod = 'KM',method = 'boot',u=36,m=80,B=1000)
plot(cal,lwd=2,lty=1,errbar.col=c(rgb(0,118,192,maxColorValue = 255)),
     xlim=c(0.0,1),ylim=c(0.0,1),
     xlab="Nomogram-predicted Probability of 3-year OS",
     ylab="Actual 5-year OS(proportion)",
     col=c(rgb(192,98,83,maxColorValue = 255)))
lines(cal[,c("mean.predicted","KM")],type = "b",lwd=2,col=c(rgb(192,98,83,maxColorValue = 255)),pch=16) 
abline(0,1,lty=3,lwd=2,col=c(rgb(0,118,192,maxColorValue = 255)))

##计算净增重新分类指数NRI
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
