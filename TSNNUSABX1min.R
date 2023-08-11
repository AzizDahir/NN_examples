setwd("C:/Users/abdulazizdahir/Desktop")
BNS<-as.matrix(read.csv("C:/Users/abdulazizdahir/Desktop/USABX1Min.csv"))
price<-as.numeric(BNS[3:nrow(BNS),2])
PV<-as.numeric(BNS[3:nrow(BNS),3])
bid<-as.numeric(BNS[3:nrow(BNS),6])
BV<-as.numeric(BNS[3:nrow(BNS),7])
ask<-as.numeric(BNS[3:nrow(BNS),10])
AV<-as.numeric(BNS[3:nrow(BNS),11])
returns<-function(aa){diff(aa)/aa[-length(aa)]}
PV.r<-returns(PV)
p.r<-returns(price)
b.r<-returns(bid)
BV.r<-returns(BV)
a.r<-returns(ask)
AV.r<-returns(AV)
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
ntr<-xx.ret[1:round(0.8*nrow(xx.ret)),]

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
pp<-VAR(ntr)$bic
pbic<-which(pp==min(pp))-1
nn<-nrow(xx.ret)
art<-NULL
for (i in 1:pbic){
art[i]<-xx.ret[-c(1:(pbic-i),nn,nn-(1:i)),3]}
art<-xx.ret[-c(1:(pbic-1),nn),3]

