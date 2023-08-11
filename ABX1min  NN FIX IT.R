#set working directory
setwd("   ")
#Read in historical stock data
BNS<-as.matrix(read.csv("  .csv"))
price<-as.numeric(BNS[3:4633,2])
bid<-as.numeric(BNS[3:4633,5])
ask<-as.numeric(BNS[3:4633,8])
DATA<-cbind(price,bid,ask)
DATA<-na.omit(DATA)



returns<-function(aa){diff(aa)/aa[-length(aa)]}
p.r<-returns(price)
b.r<-returns(bid)
a.r<-returns(ask)
xx.ret<-cbind(b.r[-1],a.r[-1],p.r[-length(p.r)],p.r[-1])
colnames(xx.ret)<-c("bidreturns","askreturns","pricemreturns","pricereturns")
xx.ret<-na.omit(xx.ret)

train<-xx.ret[1:3600,]
test<-xx.ret[3601:nrow(xx.ret),]

#maxs <- apply(xx.ret, 2, max) 
#mins <- apply(xx.ret, 2, min)

#scaled <- as.data.frame(scale(data, center = mins, scale = maxs - mins))
install.packages("neuralnet")
require(neuralnet)

n<-names(train)
f <- as.formula(paste("price returns ~", paste(n[!n %in% "price returns"], collapse = " + ")))
nn <- neuralnet(f,data=train,hidden=c(2,1),linear.output=T)


maxs <- apply(xx.ret, 2, max) 
mins <- apply(xx.ret, 2, min)

scaled <- as.data.frame(scale(xx.ret, center = mins, scale = maxs - mins))

train_ <- scaled[1:3600,]
test_ <- scaled[3601:nrow(xx.ret),]

install.packages("neuralnet")
require(neuralnet)

n <- names(train_)
f <- as.formula(paste("pricereturns ~", paste(n[!n %in% "pricereturns"], collapse = " + ")))
nn <- neuralnet(f,data=train_,hidden=c(2,1),linear.output=T)
plot(nn)
pr.nn <- compute(nn,test_[,1:3])
