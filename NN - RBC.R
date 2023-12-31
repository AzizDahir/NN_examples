# Set a seed
set.seed(500)

library(MASS)
BNS<-as.matrix(read.csv("  .csv"))
price<-as.numeric(BNS[3:4247,2])
bid<-as.numeric(BNS[3:4247,5])
ask<-as.numeric(BNS[3:4247,8])

n<-length(price)
returns<-function(aa){diff(aa)/aa[-length(aa)]}
p.r<-returns(price)
b.r<-returns(bid)
a.r<-returns(ask)
xx.ret<-cbind(b.r,p.r,a.r)

data <- xx.ret

# Check that no data is missing
apply(data,2,function(x) sum(is.na(x)))

# Train-test random splitting for linear model
#index <- sample(1:nrow(data),round(0.75*nrow(data)))
#train <- data[index,]
#test <- data[-index,]
n<-length(p.r)
p<-0.75*n
pp<-p+1
train<-xx.ret[c(1:p),]
test<-xx.ret[c(pp:n),]
train1<-data.frame(train)
test1<-data.frame(test)

# Fitting linear model
lm.fit <- glm(train1[,2]~., data=train1)
summary(lm.fit)

# Predicted data from lm
pr.lm <- predict(lm.fit,test1)

# Test MSE
MSE.lm <- sum((pr.lm - test1[,2])^2)/nrow(test1)

#-------------------------------------------------------------------------------
# Neural net fitting

# Scaling data for the NN
maxs <- apply(xx.ret, 2, max) 
mins <- apply(xx.ret, 2, min)
scaled <- as.data.frame(scale(xx.ret, center = mins, scale = maxs - mins))

# Train-test split
train_ <- scaled[index,]
test_ <- scaled[-index,]

# NN training
library(neuralnet)
n <- names(train_)
f <- as.formula(paste("train1", paste(n[!n %in% "train1"], collapse = " + ")))
nn <- neuralnet(f,data=train_,hidden=c(5,3),linear.output=T)

# Visual plot of the model
plot(nn)

# Predict
pr.nn <- compute(nn,test_[,c(1,3)])

# Results from NN are normalized (scaled)
# Descaling for comparison
xxret<-data.frame(xx.ret)
pr.nn_ <- pr.nn$net.result*(max(xxret$p.r)-min(xxret$p.r))+min(xxret$p.r)
test.r <- (test_$p.r)*(max(xxret$p.r)-min(xxret$p.r))+min(xxret$p.r)

# Calculating MSE
MSE.nn <- sum((test.r - pr.nn_)^2)/nrow(test_)

# Compare the two MSEs
print(paste(MSE.lm,MSE.nn))

# Plot predictions
par(mfrow=c(1,2))

plot(test$p.r,pr.nn_,col='red',main='Real vs predicted NN',pch=18,cex=0.7)
abline(0,1,lwd=2)
legend('bottomright',legend='NN',pch=18,col='red', bty='n')

plot(test$medv,pr.lm,col='blue',main='Real vs predicted lm',pch=18, cex=0.7)
abline(0,1,lwd=2)
legend('bottomright',legend='LM',pch=18,col='blue', bty='n', cex=.95)

# Compare predictions on the same plot
plot(test$medv,pr.nn_,col='red',main='Real vs predicted NN',pch=18,cex=0.7)
points(test$medv,pr.lm,col='blue',pch=18,cex=0.7)
abline(0,1,lwd=2)
legend('bottomright',legend=c('NN','LM'),pch=18,col=c('red','blue'))

#-------------------------------------------------------------------------------
# Cross validating

library(boot)
set.seed(200)

# Linear model cross validation
lm.fit <- glm(medv~.,data=data)
cv.glm(data,lm.fit,K=10)$delta[1]


# Neural net cross validation
set.seed(450)
cv.error <- NULL
k <- 10

# Initialize progress bar
library(plyr) 
pbar <- create_progress_bar('text')
pbar$init(k)

for(i in 1:k){
    index <- sample(1:nrow(data),round(0.9*nrow(data)))
    train.cv <- scaled[index,]
    test.cv <- scaled[-index,]
    
    nn <- neuralnet(f,data=train.cv,hidden=c(5,2),linear.output=T)
    
    pr.nn <- compute(nn,test.cv[,1:13])
    pr.nn <- pr.nn$net.result*(max(data$medv)-min(data$medv))+min(data$medv)
    
    test.cv.r <- (test.cv$medv)*(max(data$medv)-min(data$medv))+min(data$medv)
    
    cv.error[i] <- sum((test.cv.r - pr.nn)^2)/nrow(test.cv)
    
    pbar$step()
}

# Average MSE
mean(cv.error)

# MSE vector from CV
cv.error

# Visual plot of CV results
boxplot(cv.error,xlab='MSE CV',col='cyan',
        border='blue',names='CV error (MSE)',
        main='CV error (MSE) for NN',horizontal=TRUE)