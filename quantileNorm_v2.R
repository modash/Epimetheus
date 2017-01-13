suppressMessages(require(limma))
args <- commandArgs(trailingOnly = TRUE)
raw<-read.table(args[1],sep="\t",header=T,check.names=F)
nn<-ncol(raw);
#bnom<-raw[,4:nn]
anom<-normalizeQuantiles(as.matrix(raw[,4:nn]))
anom<-cbind(raw[,1:3],anom)
rm(raw)
anom<-format(anom,digits=2,trim=T,nsmall=2,scientific=F)
write.table(anom,args[2],sep="\t",row.names=F,quote=FALSE)
rm(list=ls())
