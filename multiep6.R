setwd("  ")
install.packages("neuralnet")
library(neuralnet)
install.packages("MTS")
library(MTS)
install.packages("rmgarch")
library(rmgarch)
TOTAL<-function(a,b){
BNS<-as.matrix(read.csv(" .csv"))
price<-as.numeric(BNS[3:202,2])
bid<-as.numeric(BNS[3:202,5])
ask<-as.numeric(BNS[3:202,8])
#DATA<-cbind(price,bid,ask)
#DATA<-na.omit(DATA)
#price<-as.numeric(DATA[a:b,1])
#bid<-as.numeric(DATA[a:b,2])
#ask<-as.numeric(DATA[a:b,3])
n<-length(price)
returns<-function(aa){diff(aa)/aa[-length(aa)]}
p.r<-returns(price)
b.r<-returns(bid)
a.r<-returns(ask)
xx.ret<-cbind(b.r,p.r,a.r)


VAR<-function (x, maxp = 20, output = T) 
{
    x1 = as.matrix(x)
    nT = nrow(x1)
    k = ncol(x1)
    ksq = k * k
    if (maxp < 1) 
        maxp = 1
    enob = nT - maxp
    y = x1[(maxp + 1):nT, ]
    ist = maxp + 1
    xmtx = cbind(rep(1, enob), x1[maxp:(nT - 1), ])
    if (maxp > 1) {
        for (i in 2:maxp) {
            xmtx = cbind(xmtx, x1[(ist - i):(nT - i), ])
        }
    }
    chidet = rep(0, (maxp + 1))
    s = cov(y) * (enob - 1)/enob
    chidet[1] = log(det(s))
    aic = rep(0, (maxp + 1))
    aic[1] = chidet[1]
    bic = aic
    hq = aic
    y = as.matrix(y)
    for (p in 1:maxp) {
        idm = k * p + 1
        xm = xmtx[, 1:idm]
        xm = as.matrix(xm)
        xpx <- crossprod(xm, xm)
        xpy <- crossprod(xm, y)
        beta <- solve(xpx, xpy)
        yhat <- xm %*% beta
        resi <- y - yhat
        sse <- crossprod(resi, resi)/enob
        d1 = log(det(sse))
        aic[p + 1] = d1 + (2 * p * ksq)/nT
        bic[p + 1] = d1 + (log(nT) * p * ksq)/nT
        hq[p + 1] = d1 + (2 * log(log(nT)) * p * ksq)/nT
        chidet[p + 1] = d1
    }
    maic = min(aic)
    aicor = c(1:(maxp + 1))[aic == maic] - 1
    mbic = min(bic)
    bicor = c(1:(maxp + 1))[bic == mbic] - 1
    mhq = min(hq)
    hqor = c(1:(maxp + 1))[hq == mhq] - 1
    Mstat = rep(0, maxp)
    pv = rep(0, maxp)
    for (j in 1:maxp) {
        Mstat[j] = (nT - maxp - k * j - 1.5) * (chidet[j] - chidet[j + 
            1])
        pv[j] = 1 - pchisq(Mstat[j], ksq)
    }
    if (output) {
        cat("selected order: aic = ", aicor, "\n")
        cat("selected order: bic = ", bicor, "\n")
        cat("selected order: hq = ", hqor, "\n")
    }
    if (output) {
        n1 = length(aic) - 1
        cri = cbind(c(0:n1), aic, bic, hq, c(0, Mstat), c(0, 
            pv))
        colnames(cri) <- c("p", "AIC", "BIC", "HQ", "M(p)", "p-value")
        cat("Summary table: ", "\n")
        print(round(cri, 4))
    }
    VARorder <- list(aic = aic, aicor = aicor, bic = bic, bicor = bicor, 
        hq = hq, hqor = hqor, Mstat = Mstat, Mpv = pv)
}
pp<-VAR(xx.ret)$bic
pbic<-which(pp==min(pp))-1

VMA<-function (x, lag = 20) 
{
    if (!is.matrix(x)) 
        x = as.matrix(x)
    nr = dim(x)[1]
    nc = dim(x)[2]
    g0 = var(x)
    ginv = solve(g0)
    qm = NULL
    for (i in 1:lag) {
        x1 = x[(i + 1):nr, ]
        x2 = x[1:(nr - i), ]
        g = cov(x1, x2)
        g = g * (nr - i - 1)/(nr - 1)
        h = t(g) %*% ginv %*% g %*% ginv
        qmi = nr * nr * sum(diag(h))/(nr - i)
        qm = c(qmi, qm)
    }
    tst = rev(cumsum(qm))
    ksq = nc * nc
    df = ksq * lag
    QM = NULL
    for (i in 1:lag) {
        pv = 1 - pchisq(tst[i], df)
        QM = rbind(QM, c(i, tst[i], pv))
        df = df - ksq
    }
    pvs = QM[, 3]
    dimnames(QM) = list(names(pvs), c("  j  ", "  Q(j,m) ", " p-value"))
    cat("Q(j,m) Statistics: ", "\n")
    print(QM[,3])
    #printCoefmat(QM, digits = 3)
    #par(mfcol = c(1, 1))
    #plot(pvs, ylim = c(0, 1), xlab = "j", ylab = "prob", main = "p-values: Q(j,m) Statistics")
    #abline(h = c(0))
    #lines(rep(0.05, lag), lty = 2, col = "blue")
}
qq<-NULL
for (i in 1:20){
qq[i]<-tail(VMA(xx.ret,i),1)}
#qq<-VMA(xx.ret,4)
qq1<-max.col(t(qq) > 0.05, "first")
q<-qq1-1

garch11.spec = ugarchspec(mean.model = list(armaOrder = c(pbic,q)),
variance.model = list(garchOrder = c(1,1),
                          model = "sGARCH"),
                          distribution.model = "norm")

dcc.garch11.spec = dccspec(uspec = multispec(replicate(3, garch11.spec) ),
                           dccOrder = c(1,1),
                           distribution = "mvnorm")
dcc.garch11.spec

dcc.fit = dccfit(dcc.garch11.spec, data = xx.ret )
n<-length(p.r)
err<-function(i){dcc.fit@model$residuals[i,2]}
errb<-function(i){dcc.fit@model$residuals[i,1]}
erra<-function(i){dcc.fit@model$residuals[i,3]}
berr<-errb(1:n)
perr<-err(1:n)
aerr<-erra(1:n)
system.time(fitv<-VARMACpp(xx.ret, p = pbic, q = q, include.mean = T,
fixed = NULL, beta=NULL, sebeta=NULL,
prelim = F, details = F, thres = 2))

#Problems
#g<-fitv$Phi[2,]
g<-matrix(fitv$Phi[2,],nrow=pbic,byrow=T)
bbp<-g[,2]
bbb<-g[,1]
bba<-g[,3]
#bb<-g[seq(2, length(g), 3)]
#bb<-matrix(g,ncol=3,byrow=T)

h<-fitv$Theta[2,]
#h<-matrix(fitv$Theta[2,],nrow=q,byrow=T)
#aa<-function(i){g[i,2]}
#aa<-h[,2]
aap<-h[seq(2, length(h), 3)]
aab<-h[seq(1, length(h), 3)]
aaa<-h[seq(3, length(h), 3)]
#aa<-matrix(h,ncol=3,byrow=T)

####Building the model
AR<-NULL
for (i in 1:pbic){
AR[i]<-bbp[i]*tail(p.r,i)[1]+bbb[i]*tail(b.r,i)[1]+bba[i]*tail(a.r,i)[1]}
arr<-sum(AR)
MA<-NULL
for (i in 1:q){
#MA[i]<-aap[i]*tail(rrep,i)[1]+aab[i]*tail(rreb,i)[1]+aaa[i]*tail(rrea,i)[1]}
MA[i]<-aap[i]*tail(resids[,2],i)[1]+aab[i]*tail(resids[,1],i)[1]+aaa[i]*tail(resids[,3],i)[1]}

maa<-sum(MA)
rtt<-arr+maa
ep<-(1+rtt)*tail(price,1)[1]
return(ep)}
predp0<-p.r[3:n]
predp1<-bbp[1]*p.r[2:c(n-1)]
predp2<-bbp[2]*p.r[1:c(n-2)]
per<-aap*perr[2:c(n-1)]
predb1<-bbb[1]*b.r[2:c(n-1)]
predb2<-bbb[2]*b.r[1:c(n-2)]
ber<-aab*berr[2:c(n-1)]
preda1<-bba[1]*a.r[2:c(n-1)]
preda2<-bba[2]*a.r[1:c(n-2)]
aer<-aaa*aerr[2:c(n-1)]
xrt<-cbind(predp0, predp1,predp2,per,predb1,predb2,ber,preda1,preda2,aer)
f<-predp0~ predp1 + predp2 + per + predb1 + predb2 + ber + preda1 + preda2 + aer
#f <- as.formula(paste("medv ~", paste(n[!n %in% "medv"], collapse = " + ")))
nn <- neuralnet(f,data=xrt,hidden=c(5,3),linear.output=T)
plot(nn)
pr.nn <- compute(nn,xrt[,2:10])
xxret<-data.frame(xrt)
pr.nn_ <- pr.nn$net.result*(max(xxret$predp0)-min(xxret$predp0))+min(xxret$predp0)
#test.r <- (xrt$predp0)*(max(xxret$predp0)-min(xxret$predp0))+min(xxret$predp0)

nnep<-(1+tail(pr.nn_,1))*tail(price,1)[1]