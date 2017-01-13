args <- commandArgs(trailingOnly = TRUE)
print(args)
png(args[3])
title<-paste(args[4]," - MA Plot:",args[5],"_Vs_",args[6], " - ", args[7])
denom<-read.table(args[1],sep="\t",header=F,skip=1)
numer<-read.table(args[2],sep="\t",header=F,skip=1)
denom<-denom[,4]
numer<-numer[,4]
denom[denom<=1]<-1
numer[numer<=1]<-1
M<-log2(numer/denom)
A<-0.5*log2(numer*denom)
M <- as.matrix(M)
A <- as.matrix(A)
fitbrut<-lowess(A,M,f=2/3)
plot(A,M,main=title,col=rgb(100,0,0,80,maxColorValue=255),pch=20)
abline(h=0,col="black",lty=5)
lines(fitbrut,col="green")
dev.off()
rm(list=ls())
