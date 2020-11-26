## histogramme.R
## 

# input
cmd=commandArgs(TRUE)
annotFile=cmd[1]
file=cmd[2]
name=cmd[3]
lgMin=as.integer(cmd[4])
lgMax=as.integer(cmd[5])

#annot
dat   <- read.table(annotFile, skip = 2, header = FALSE, sep ='\t')
myAnnot=as.vector(dat$V1)
myAnnot=append(myAnnot,"non-annotated")

#
colFil="Araport11_Annot/annot_colors.txt"
coldat=read.table(colFil,header=T,sep="\t")
myCol=NULL
for (i in 1:length(myAnnot)){
tmpAnnot=myAnnot[i]
myCol=append(myCol,as.vector(coldat[coldat$annot==tmpAnnot,]$color))
}
#
data=read.table(file,header=T)
dataFilter=data[data$lg >=lgMin & data$lg <= lgMax,]

#mat=as.matrix(data[6:26,2:13])
mat=as.matrix(dataFilter[,1:length(myAnnot)+1])
# output
pdf(file=paste(name,".pdf",sep=""),width=8,height=8)
layout(matrix(c(1,1,2),ncol=3))
par(adj=0.5,bty="n",las=0,lab=c(20,10,10))
barplot(t(mat),col=myCol,names.arg=seq(lgMin,lgMax,1),main=name,ylab="Number of reads",xlab="Reads length")

tmp=matrix(0,ncol=length(colnames(mat)))
barplot(tmp,col=myCol,legend=c(colnames(mat)),yaxt="n",ylab="",xlab="",axes=FALSE,border='NA')
invisible(dev.off())
