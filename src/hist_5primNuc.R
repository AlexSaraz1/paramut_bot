

cmd=commandArgs(TRUE)
infile=cmd[1]
outname=cmd[2]
lgMin=as.integer(cmd[3])
lgMax=as.integer(cmd[4])


myCol=c("dark green","dark red","darkgoldenrod1","dark blue","grey")		

data=read.table(infile,header=T,sep="\t")


pdf(file=paste(outname,".pdf",sep=""),width=8,height=8)
barplot(t(as.matrix(data)),col=myCol,legend=colnames(data),names.arg=seq(lgMin,lgMax,1),ylab="Number of reads",xlab="Reads length")
invisible(dev.off())

databis=matrix(0,nrow=length(data[,1]),ncol=length(data[1,]))
for (r in 1:length(data[,1])){
databis[r,]=as.vector(unlist(data[r,]*100/sum(data[r,])))
}

mat=as.matrix(data)
pdf(file=paste("proportion_",outname,".pdf",sep=""),width=8,height=8)
layout(matrix(c(1,1,2),ncol=3))
par(adj=0.5,bty="n",las=0,lab=c(20,10,10))
barplot(t(as.matrix(databis)),col=myCol,names.arg=seq(lgMin,lgMax,1),ylab="Proportion of reads",xlab="Reads length")
tmp=matrix(0,ncol=length(colnames(mat)))
barplot(t(tmp),col=myCol,yaxt="n",ylab="",xlab="",axes=FALSE,border='NA')
legend("topleft",legend =c(colnames(mat)),cex=2,fill=myCol)
invisible(dev.off())


