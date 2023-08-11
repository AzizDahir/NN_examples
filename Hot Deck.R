setwd(" ")
require(mlbench)
require(car)
require(MASS)
require(plyr)
data("PimaIndiansDiabetes")
PID<-PimaIndiansDiabetes
id<-1:nrow(PID)
PID[,9]<-as.numeric(PID[,9])
## 2 -> pos 1 -> neg
PID.rs<-rowSums(PID)
#n.PID<-cbind(PID,PID.rs,id)
n.PID<-cbind(PID,PID.rs)
n.PID$glucose[n.PID$glucose==0]<-NA
n.PID[,2]<-as.numeric(n.PID[,2])
n.PID$pressure[n.PID$pressure==0]<-NA
n.PID[,3]<-as.numeric(n.PID[,3])
n.PID$triceps[n.PID$triceps==0]<-NA
n.PID[,4]<-as.numeric(n.PID[,4])
n.PID$insulin[n.PID$insulin==0]<-NA
n.PID[,5]<-as.numeric(n.PID[,5])
n.PID$mass[n.PID$mass==0]<-NA
n.PID[,6]<-as.numeric(n.PID[,6])
order.PID<-n.PID[order(n.PID[,10]),]
PICK<-function(x,y){order.PID[c((which(order.PID$id==x)-y):(which(order.PID$id==x)+y)),]}
PICK(768,20)

#Hot Deck imputation
order.PID1<-order.PID
for (i in 1:ncol(order.PID)){
  LIS<-which(is.na(order.PID[,i]))
  if(length(LIS)!=0){
    for (j in 1:length(LIS)){
      new.PID<-cbind(order.PID,abs(order.PID$PID.rs[LIS[j]]-order.PID$PID.rs))
      new.PID<-new.PID[order(abs(order.PID$PID.rs[LIS[j]] - order.PID$PID.rs)),]
        for (k in 2:nrow(PID)){if(!is.na(new.PID[k,i])){new.PID[1,i]<-new.PID[k,i];break}}
      order.PID1[LIS[j],i]<-new.PID[1,i]
    }}}
View(order.PID1)


### Clustering + Logistic regression to see how good it is to see to predict if missing data
C.PID<-na.omit(n.PID)
cor(C.PID)

C.PID0<-subset(C.PID,C.PID$diabetes==1)
summary(C.PID0)
C.PID1<-subset(C.PID,C.PID$diabetes==2)
summary(C.PID1)

flag.glu<-is.na(n.PID$glucose)
flag.glu<-as.numeric(flag.glu)
flag.pre<-is.na(n.PID$pressure)
flag.pre<-as.numeric(flag.pre)
flag.tri<-is.na(n.PID$triceps)
flag.tri<-as.numeric(flag.tri)
flag.ins<-is.na(n.PID$insulin)
flag.ins<-as.numeric(flag.ins)
flag.mass<-is.na(n.PID$mass)
flag.mass<-as.numeric(flag.mass)
# 1 missing 0 not missing
n.PID<-cbind(n.PID,flag.glu,flag.pre,flag.tri,flag.ins,flag.mass)
n.PID<-n.PID[c(1,2,11,3,12,4,13,5,14,6,15,7:10)]

diabetes<-n.PID$diabetes-1
n.Pid<-cbind(n.PID,diabetes)
n.PID<-n.Pid[,-14]

n.PID$flag.glu[n.PID$flag.glu==0]<-5
n.PID$flag.glu[n.PID$flag.glu==1]<-4
flag.gluc<-n.PID$flag.glu-4
n.PID<-n.PID[,-3]
n.PID<-cbind(n.PID,flag.gluc)
n.PID<-n.PID[c(1:2,15,3:14)]

n.PID$flag.pre[n.PID$flag.pre==0]<-5
n.PID$flag.pre[n.PID$flag.pre==1]<-4
flag.pres<-n.PID$flag.pre-4
n.PID<-n.PID[,-5]
n.PID<-cbind(n.PID,flag.pres)
n.PID<-n.PID[c(1:4,15,5:14)]

n.PID$flag.tri[n.PID$flag.tri==0]<-5
n.PID$flag.tri[n.PID$flag.tri==1]<-4
flag.tric<-n.PID$flag.tri-4
n.PID<-n.PID[,-7]
n.PID<-cbind(n.PID,flag.tric)
n.PID<-n.PID[c(1:6,15,7:14)]

n.PID$flag.ins[n.PID$flag.ins==0]<-5
n.PID$flag.ins[n.PID$flag.ins==1]<-4
flag.insu<-n.PID$flag.ins-4
n.PID<-n.PID[,-9]
n.PID<-cbind(n.PID,flag.insu)
n.PID<-n.PID[c(1:8,15,9:14)]

n.PID$flag.mass[n.PID$flag.mass==0]<-5
n.PID$flag.mass[n.PID$flag.mass==1]<-4
flag.Mass<-n.PID$flag.mass-4
n.PID<-n.PID[,-11]
n.PID<-cbind(n.PID,flag.Mass)
n.PID<-n.PID[c(1:10,15,11:14)]

indexes<-sample(1:nrow(n.PID),size=ceiling(0.4*nrow(n.PID)))
test<-n.PID[indexes,]
training<-n.PID[-indexes,]

##glu.glm<-glm(is.na(glucose)~pressure+triceps+pregnant+mass+diabetes+
###pedigree+age,data=training,family=binomial(link="logit"))
##summary(glu.glm)
### nothing, nothing here is predictive to why data is missing
##PRE.glu<-predict(glu.glm,test,type="response")
##(a.glu<-table(PRE.glu>0.5,is.na(test$glucose)))

glu.glm<-glm(is.na(glucose)~diabetes
,data=training,family=binomial(link="logit"))
summary(glu.glm)
### nothing, nothing here is predictive to why data is missing
PRE.glu<-predict(glu.glm,test,type="response")
(a.glu<-table(PRE.glu>0.5,is.na(test$glucose)))


C.Pid.glu<-lm(glucose~diabetes+insulin+triceps+pedigree+
mass+pregnant+age,data=C.Pid)
summary(C.Pid.glu)
vif(C.Pid.glu)
step <- stepAIC(C.Pid.glu, direction="both")
step$anova


pre.glm<-glm(is.na(pressure)~mass+pregnant+glucose+triceps+diabetes+
pedigree+age,data=training,family=binomial)
summary(pre.glm)
#age mass
PRE.pre<-predict(pre.glm,test,type="response")
a<-table(PRE.pre>0.5,is.na(test$pressure))

tri.glm<-glm(is.na(triceps)~	mass+
				      pressure+
					pregnant+
					glucose+
					diabetes+
					pedigree+
					age,
data=n.PID,family=binomial(link="logit"))
summary(tri.glm)
### mass pressure pedigree age
PRE<-predict(tri.glm,test,type="response")
a<-table(PRE>0.5,is.na(test$triceps))
(a[1,1]+a[2,2])/sum(rowSums(a))
## 74% accurate
##Recall for false (not missing) 
a[1,1]/rowSums(a)[1]  
##Precision for false column
a[2,1]/colSums(a)[1]
#Recall for true (missing)
a[2,2]/rowSums(a)[2]
##Precision for true column
a[1,2]/colSums(a)[2]

ins.glm<-glm(is.na(insulin)~mass+glucose+pressure+triceps+pregnant+diabetes+
pedigree+age,data=training,family=binomial(link="logit"))
summary(ins.glm)
### age
PRE.ins<-predict(ins.glm,test,type="response")
(a.ins<-table(PRE.ins>0.5,is.na(test$insulin)))
(a.ins[1,1]+a.ins[2,2])/sum(rowSums(a.ins))
## 71% accuracy
##Recall for false (not missing) 
a.ins[1,1]/rowSums(a.ins)[1] 
##Precision for false column
a.ins[2,1]/colSums(a.ins)[1]
#Recall for true (missing)
a.ins[2,2]/rowSums(a.ins)[2]
##Precision for true column
a.ins[1,2]/colSums(a.ins)[2]

mass.glm<-glm(is.na(mass)~glucose+pregnant+diabetes+pedigree+pressure+
age,data=n.PID,family=binomial(link="logit"))
summary(mass.glm)
#nothing maybe pressure
PRE.mass<-predict(mass.glm,test,type="response")
(a.mass<-table(PRE.mass>0.5,is.na(test$mass)))

####################################
Casewise Deletion
######################################
C.PID<-na.omit(n.PID)
C.Pid<-C.PID[,c(-3,-5,-7,-9,-11)]
summary(N.pid,use="pairwise.complete.obs") 
## insulin,glucose,row sums high correlation
#glucose and diabetes
## triceps mass
C.PID0<-subset(C.Pid,C.Pid$diabetes==0)
summary(C.PID0)
C.PID1<-subset(C.Pid,C.PID$diabetes==1)
summary(C.PID1)
########################################################
Pairwise
#########################################################
N.pid<-n.PID[,c(-3,-5,-7,-9,-11,-14)]
cor(N.pid,use="pairwise.complete.obs") 
summary(N.pid,use="pairwise.complete.obs") 
mean(N.pid,na.rm=TRUE)
#########################################################
Mean imputation
#########################################################
id<-1:nrow(n.PID)
n.PID<-cbind(n.PID,id)
impute.mean <- function(x) replace(x, is.na(x), mean(x, na.rm = TRUE))
n.PID2 <- ddply(n.PID, ~diabetes, transform, glucose = impute.mean(glucose),
     pressure = impute.mean(pressure),triceps=impute.mean(triceps),
insulin=impute.mean(insulin),mass=impute.mean(mass))
n.PID2<-n.PID2[order(n.PID2$id), ]
###diabetic and non diabetic have different means
#################################################################
Medium Imputations
#################################################################
impute <- function(x, fun) {
 missing <- is.na(x)
  replace(x, missing, fun(x[!missing]))
}
n.PID3<-ddply(n.PID, ~ diabetes, transform, glucose = impute(glucose,median),
     pressure = impute(pressure,median),triceps=impute(triceps,median),
insulin=impute(insulin,median),mass=impute(mass,median))
n.PID3<-n.PID3[order(n.PID3$id), ]
##################################################################
Min Imputation
##################################################################
n.PID4<-ddply(n.PID, ~ diabetes, transform, glucose = impute(glucose,min),
     pressure = impute(pressure,min),triceps=impute(triceps,min),
insulin=impute(insulin,min),mass=impute(mass,min))
n.PID4<-n.PID4[order(n.PID4$id), ]
##################################################################
Max Imputation
##################################################################
n.PID5<-ddply(n.PID, ~ diabetes, transform, glucose = impute(glucose,max),
     pressure = impute(pressure,max),triceps=impute(triceps,max),
insulin=impute(insulin,max),mass=impute(mass,max))
n.PID5<-n.PID5[order(n.PID5$id), ]
###################################################################

C.Pid.glm<-glm(diabetes~pedigree+glucose+mass+pregnant+
age,data=C.Pid,family=binomial)
summary(C.Pid.glm)
vif(C.Pid.glm)
step <- stepAIC(C.Pid.glm, direction="both")
step$anova
######################################################
#########################################################
pre.lm<-lm(is.na(pressure)~diabetes+age+mass+pregnant,data=training)
vif(pre.lm)

pre.glm<-glm(is.na(pressure)~diabetes+age+mass+pregnant
,data=training,family=binomial)
summary(pre.glm)
PRE.pre<-predict(pre.glm,test,type="response")
(a.pre<-table(PRE.pre>0.5,is.na(test$pressure)))
vif(pre.glm)

pre.glm<-glm(is.na(pressure)~pregnant+pedigree+
diabetes,data=training,family=binomial)
summary(pre.glm)

n.PID0<-subset(n.PID,n.PID$diabetes==0)
n.PID1<-subset(n.PID,n.PID$diabetes==1)


age1.PID0<-subset(n.PID0,n.PID0$age<24)
age2.PID0<-subset(n.PID0,n.PID0$age>23 & age<28)
age3.PID0<-subset(n.PID0,n.PID0$age>27 & age<38)
age4.PID0<-subset(n.PID0,n.PID0$age>37)

order0.age1<-age1.PID0[order(age1.PID0[,14]),]
order0.age2<-age2.PID0[order(age2.PID0[,14]),]
order0.age3<-age3.PID0[order(age3.PID0[,14]),]
order0.age4<-age4.PID0[order(age4.PID0[,14]),]

age1.PID1<-subset(n.PID1,n.PID1$age<29)
age2.PID1<-subset(n.PID1,n.PID1$age>28 & age<37)
age3.PID1<-subset(n.PID1,n.PID1$age>36 & age<45)
age4.PID1<-subset(n.PID1,n.PID1$age>44)

order1.age1<-age1.PID1[order(age1.PID1[,14]),]
order1.age2<-age2.PID1[order(age2.PID1[,14]),]
order1.age3<-age3.PID1[order(age3.PID1[,14]),]
order1.age4<-age4.PID1[order(age4.PID1[,14]),]

indexes<-sample(1:nrow(n.PID0),size=ceiling(0.4*nrow(n.PID0)))
test0<-n.PID0[indexes,]
training0<-n.PID0[-indexes,]

ins0.glm<-glm(is.na(insulin)~mass+pressure+triceps+pregnant+
pedigree+age+glucose,data=training0,family=binomial(link="logit"))
summary(ins0.glm)
vif(ins0.glm) 
PRE0.ins<-predict(ins0.glm,test,type="response")
(a.ins0<-table(PRE0.ins>0.5,is.na(test$insulin)))
####multicollinearity preg and age
### got rid of diabetes sinces we classified for diabetes (always 0)
ins0.glm<-glm(is.na(insulin)~pressure+triceps+
pedigree+glucose,data=training0,family=binomial(link="logit"))
summary(ins0.glm)
vif(ins0.glm) 
PRE0.ins<-predict(ins0.glm,test0,type="response")
(a.ins0<-table(PRE0.ins>0.5,is.na(test0$insulin)))

(a.ins0[1,1]+a.ins0[2,2])/sum(rowSums(a.ins0))
## 72% accuracy
##Recall for false (not missing) 
a.ins0[1,1]/rowSums(a.ins0)[1] 
##71% correctly classified as false (not missing)
##Precision for false column
a.ins0[2,1]/colSums(a.ins0)[1]
## 0% classified as incorrectly as false
#Recall for true (missing)
a.ins0[2,2]/rowSums(a.ins0)[2]
## all of the predicted as missing are 100%
##Precision for true column
a.ins0[1,2]/colSums(a.ins0)[2]
##95% incorrectly labelled as missing when truly not missing 
########################################################################

tri0.glm<-glm(is.na(triceps)~	mass+
				      pressure+
					pregnant+
					glucose+
					pedigree+
					age,
data=training0,family=binomial(link="logit"))
summary(tri.glm)
vif(tri0.glm)
### significant -> mass pressure pedigree age
PRE.tri0<-predict(tri0.glm,test0,type="response")
(a.tri0<-table(PRE.tri0>0.5,is.na(test0$triceps)))

(a.tri0[1,1]+a.tri0[2,2])/sum(rowSums(a.tri0))
## 74% accuracy
##Recall for false (not missing) 
a.tri0[1,1]/rowSums(a.tri0)[1] 
##78% correctly classified as false (not missing)
##Precision for false column
a.tri0[2,1]/colSums(a.tri0)[1]
## 10% classified as incorrectly as false
#Recall for true (missing)
a.tri0[2,2]/rowSums(a.tri0)[2]
## all of the predicted as missing are 42%
##Precision for true column
a.ins0[1,2]/colSums(a.ins0)[2]
##97% incorrectly labelled as missing when truly not missing 
##########################################################################

mass0.glm<-glm(is.na(mass)~glucose+pregnant+pedigree+pressure+
age,data=training0,family=binomial(link="logit"))
step <- stepAIC(mass0.glm, direction="both")
step$anova
summary(mass0.glm)
vif(mass0.glm)

PRE.mass0<-predict(mass.glm,test0,type="response")
(a.mass0<-table(PRE.mass0>0.025,is.na(test0$mass)))

summary(PRE.mass0)
##########################################################################
order0.age1<-age1.PID0[order(age1.PID0[,14]),]
order0.age2<-age2.PID0[order(age2.PID0[,14]),]
order0.age3<-age3.PID0[order(age3.PID0[,14]),]
order0.age4<-age4.PID0[order(age4.PID0[,14]),]

order1.age1<-age1.PID1[order(age1.PID1[,14]),]
order1.age2<-age2.PID1[order(age2.PID1[,14]),]
order1.age3<-age3.PID1[order(age3.PID1[,14]),]
order1.age4<-age4.PID1[order(age4.PID1[,14]),]

##########################################################################
Hot Deck Imputation
##########################################################################
Order0.PID2<-order0.age2
for (i in 1:ncol(order0.age2)){
  LIS<-which(is.na(order0.age2[,i]))
  if(length(LIS)!=0){
    for (j in 1:length(LIS)){
      new.PID<-cbind(order0.age2,abs(order0.age2$PID.rs[LIS[j]]-order0.age2$PID.rs))
      new.PID<-new.PID[order(abs(order0.age2$PID.rs[LIS[j]] - order0.age2$PID.rs)),]
        for (k in 2:nrow(PID)){if(!is.na(new.PID[k,i])){new.PID[1,i]<-new.PID[k,i];break}}
      Order0.PID2[LIS[j],i]<-new.PID[1,i]
    }}}
View(Order0.PID2)
############################################################################
Order0.PID1<-order0.age1
for (i in 1:ncol(order0.age1)){
  LIS<-which(is.na(order0.age1[,i]))
  if(length(LIS)!=0){
    for (j in 1:length(LIS)){
      new.PID<-cbind(order0.age1,abs(order0.age1$PID.rs[LIS[j]]-order0.age1$PID.rs))
      new.PID<-new.PID[order(abs(order0.age1$PID.rs[LIS[j]] - order0.age1$PID.rs)),]
        for (k in 2:nrow(PID)){if(!is.na(new.PID[k,i])){new.PID[1,i]<-new.PID[k,i];break}}
      Order0.PID1[LIS[j],i]<-new.PID[1,i]
    }}}
View(Order0.PID1)
##############################################################################
Order0.PID3<-order0.age3
for (i in 1:ncol(order0.age3)){
  LIS<-which(is.na(order0.age3[,i]))
  if(length(LIS)!=0){
    for (j in 1:length(LIS)){
      new.PID<-cbind(order0.age3,abs(order0.age3$PID.rs[LIS[j]]-order0.age3$PID.rs))
      new.PID<-new.PID[order(abs(order0.age3$PID.rs[LIS[j]] - order0.age3$PID.rs)),]
        for (k in 2:nrow(PID)){if(!is.na(new.PID[k,i])){new.PID[1,i]<-new.PID[k,i];break}}
      Order0.PID3[LIS[j],i]<-new.PID[1,i]
    }}}
View(Order0.PID3)
###############################################################################
Order0.PID4<-order0.age4
for (i in 1:ncol(order0.age4)){
  LIS<-which(is.na(order0.age4[,i]))
  if(length(LIS)!=0){
    for (j in 1:length(LIS)){
      new.PID<-cbind(order0.age4,abs(order0.age4$PID.rs[LIS[j]]-order0.age4$PID.rs))
      new.PID<-new.PID[order(abs(order0.age4$PID.rs[LIS[j]] - order0.age4$PID.rs)),]
        for (k in 2:nrow(PID)){if(!is.na(new.PID[k,i])){new.PID[1,i]<-new.PID[k,i];break}}
      Order0.PID4[LIS[j],i]<-new.PID[1,i]
    }}}
View(Order0.PID4)
#################################################################################
Order1.PID1<-order1.age1
for (i in 1:ncol(order1.age1)){
  LIS<-which(is.na(order1.age1[,i]))
  if(length(LIS)!=0){
    for (j in 1:length(LIS)){
      new.PID<-cbind(order1.age1,abs(order1.age1$PID.rs[LIS[j]]-order1.age1$PID.rs))
      new.PID<-new.PID[order(abs(order1.age1$PID.rs[LIS[j]] - order1.age1$PID.rs)),]
        for (k in 2:nrow(PID)){if(!is.na(new.PID[k,i])){new.PID[1,i]<-new.PID[k,i];break}}
      Order1.PID1[LIS[j],i]<-new.PID[1,i]
    }}}
View(Order1.PID1)
################################################################################
Order1.PID2<-order1.age2
for (i in 1:ncol(order1.age2)){
  LIS<-which(is.na(order1.age2[,i]))
  if(length(LIS)!=0){
    for (j in 1:length(LIS)){
      new.PID<-cbind(order1.age2,abs(order1.age2$PID.rs[LIS[j]]-order1.age2$PID.rs))
      new.PID<-new.PID[order(abs(order1.age2$PID.rs[LIS[j]] - order1.age2$PID.rs)),]
        for (k in 2:nrow(PID)){if(!is.na(new.PID[k,i])){new.PID[1,i]<-new.PID[k,i];break}}
      Order1.PID2[LIS[j],i]<-new.PID[1,i]
    }}}
View(Order1.PID2)
#################################################################################
Order1.PID3<-order1.age3
for (i in 1:ncol(order1.age3)){
  LIS<-which(is.na(order1.age3[,i]))
  if(length(LIS)!=0){
    for (j in 1:length(LIS)){
      new.PID<-cbind(order1.age3,abs(order1.age3$PID.rs[LIS[j]]-order1.age3$PID.rs))
      new.PID<-new.PID[order(abs(order1.age3$PID.rs[LIS[j]] - order1.age3$PID.rs)),]
        for (k in 2:nrow(PID)){if(!is.na(new.PID[k,i])){new.PID[1,i]<-new.PID[k,i];break}}
      Order1.PID3[LIS[j],i]<-new.PID[1,i]
    }}}
View(Order1.PID3)
################################################################################
Order1.PID4<-order1.age4
for (i in 1:ncol(order1.age4)){
  LIS<-which(is.na(order1.age4[,i]))
  if(length(LIS)!=0){
    for (j in 1:length(LIS)){
      new.PID<-cbind(order1.age4,abs(order1.age4$PID.rs[LIS[j]]-order1.age4$PID.rs))
      new.PID<-new.PID[order(abs(order1.age4$PID.rs[LIS[j]] - order1.age4$PID.rs)),]
        for (k in 2:nrow(PID)){if(!is.na(new.PID[k,i])){new.PID[1,i]<-new.PID[k,i];break}}
      Order1.PID4[LIS[j],i]<-new.PID[1,i]
    }}}
View(Order1.PID4)
################################################################################


