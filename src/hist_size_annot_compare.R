cmd=commandArgs(TRUE)
annotFile=cmd[1]
file=cmd[2]
name=cmd[3]
lgMin=as.integer(cmd[4])
lgMax=as.integer(cmd[5])

AnnotOfInterest=cmd[6]


colFil="Araport11_Annot/annot_colors.txt"
coldat=read.table(colFil,header=T,sep="\t")
myCol=as.vector(coldat[coldat$annot==AnnotOfInterest,]$color)

dat   <- read.table(annotFile, skip = 2, header = FALSE, sep ='\t')
myAnnot=as.vector(dat$V1)
myAnnot=append(myAnnot,"non-annotated")

data=read.table(file,header=T)
dataFilter=data[data$lg >=lgMin & data$lg <= lgMax,]


mat=matrix(0,ncol=length(seq(lgMin,lgMax,1)),nrow=2)

for (i in 1:length(seq(lgMin,lgMax,1))){
mat[1,i]=sum(dataFilter[i,c(1:length(myAnnot)+1)])
mat[2,i]=dataFilter[[AnnotOfInterest]][i]
}


mat2=matrix(0,ncol=length(seq(lgMin,lgMax,1)),nrow=2)
for (i in 1:length(seq(lgMin,lgMax,1))){
mat2[1,i]=mat[1,i]*100/sum(mat[1,])
mat2[2,i]=mat[2,i]*100/sum(mat[2,])
}

pdf(file=paste(name,"_compare",".pdf",sep=""),width=8,height=8)
par(adj=0.5,bty="n",las=0,lab=c(20,10,10))
barplot(mat2,beside=T,ylab="Proportion of reads",xlab="Reads length",legend=c("All",AnnotOfInterest),col=c("dark grey",myCol),names.arg=seq(lgMin,lgMax,1),space=c(0,0.5))
invisible(dev.off())
