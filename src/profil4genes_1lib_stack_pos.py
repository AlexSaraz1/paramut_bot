#!/usr/bin/python


import sys,os,string


geneFileId=sys.argv[1]

geneFilegff=sys.argv[2]
#
sRNAfileGff1=sys.argv[3]
libName1=sys.argv[4]
nbReads1=sys.argv[5]

# read genes id
fich=open(geneFileId,'r')
lign=fich.readline()
myGenesId={}
while lign != "" :
	lignSplit=string.split(lign.replace("\n","").replace('"',''),sep="\t")
	myId=lignSplit[0]
	myGenesId[myId]=[myId]
	print myId
	lign=fich.readline()

fich.close()


fich=open(geneFilegff,'r')
lign=fich.readline()
while lign != "" :
	lignSplit=string.split(lign.replace("\n",""),sep="\t")
	attrib=lignSplit[8]
	myId=string.split(attrib.replace("ID=",""),sep=";")[0]
	if myGenesId.get(myId) :
		#print myId
		chro=lignSplit[0]
		st=lignSplit[3]
		sp=lignSplit[4]
		## reads lib1
		awkfich=open(myId+"_cmd.awk",'w')
		awkfich.write("/^"+chro+"/{if($4 >= "+st+" && $5 <= "+sp+"){ print $0}}\n")
		awkfich.close()
		os.system("cat "+sRNAfileGff1+" | awk -f "+myId+"_cmd.awk  > reads"+libName1+"_"+myId+".gff")
		# profil file lib1
		os.system("python ./src/sRNAdistri_plusStrand_gff3.py reads"+libName1+"_"+myId+".gff " +myId+" "+st+" "+sp+" > nbReadsStranded"+libName1+"_"+myId+".txt")
		# read length lib1
		awkfich=open(myId+"_cmd.awk",'w')
		awkfich.write("/^/{ cmp=0; while (cmp < $13){ print length($NF)\"\\t\"$16;cmp=cmp+1}}END{ print 35\"\\t\"1}\n")
		awkfich.close()
		os.system("cat reads"+libName1+"_"+myId+".gff | awk -f "+myId+"_cmd.awk  > ReadLength_"+libName1+"_"+myId+".txt")
		#
		Rfich=open(myId+"_cmd.R",'w')
		# matrix for read length lib1
		Rfich.write("dt=read.table(\"ReadLength_"+libName1+"_"+myId+".txt\",header=F,sep=\"\t\")"+"\n")
		Rfich.write("dtUniq=dt[dt$V2 == 1,]"+"\n")
		Rfich.write("dtMultipl=dt[dt$V2 != 1,]"+"\n")
		Rfich.write("nbReads="+nbReads1+"\n")
		Rfich.write("countUniq=hist(dtUniq$V1,breaks=seq(14.5,35.5,1),plot = FALSE)$count*10000000/nbReads"+"\n")
		Rfich.write("countMultipl=hist(dtMultipl$V1,breaks=seq(14.5,35.5,1),plot = FALSE)$count*10000000/nbReads"+"\n")
		Rfich.write("matNbReads1=matrix(c(countUniq,countMultipl),nrow=2,byrow=T)"+"\n")
		# data1
		Rfich.write("data=read.table(\"nbReadsStranded"+libName1+"_"+myId+".txt\",header=F,sep=\"\t\")"+"\n")
		Rfich.write("nbReads="+nbReads1+"\n")
		Rfich.write("myVect21plus=data[,1]*10000000/nbReads"+"\n")
		Rfich.write("myVect22plus=data[,2]*10000000/nbReads"+"\n")
		Rfich.write("myVect24plus=data[,3]*10000000/nbReads"+"\n")
		Rfich.write("myVectOtherplus=data[,4]*10000000/nbReads"+"\n")
		Rfich.write("myVect21minus=-1*data[,5]*10000000/nbReads"+"\n")
		Rfich.write("myVect22minus=-1*data[,6]*10000000/nbReads"+"\n")
		Rfich.write("myVect24minus=-1*data[,7]*10000000/nbReads"+"\n")
		Rfich.write("myVectOtherminus=-1*data[,8]*10000000/nbReads"+"\n")
		Rfich.write("mycol=c(\"grey\",\"blue1\",\"chartreuse\",\"firebrick1\")"+"\n")
		Rfich.write("mat1=matrix(0,ncol=length(myVect21plus),nrow=8)"+"\n")
		Rfich.write("mat1[1,]=myVectOtherplus"+"\n")
		Rfich.write("mat1[2,]=myVect21plus"+"\n")
		Rfich.write("mat1[3,]=myVect22plus"+"\n")
		Rfich.write("mat1[4,]=myVect24plus"+"\n")
		Rfich.write("mat2=matrix(0,ncol=length(myVect21plus),nrow=8)"+"\n")
		Rfich.write("mat2[1,]=myVectOtherminus"+"\n")
		Rfich.write("mat2[2,]=myVect21minus"+"\n")
		Rfich.write("mat2[3,]=myVect22minus"+"\n")
		Rfich.write("mat2[4,]=myVect24minus"+"\n")
		Rfich.write("mx=max(c(myVect21plus+myVect22plus+myVect24plus+myVectOtherplus))"+"\n")
		Rfich.write("mn=min(c(myVect21minus+myVect22minus+myVect24minus+myVectOtherminus))"+"\n")
		# start graph
		Rfich.write("pdf(file=\"distriSur"+myId+".pdf\")"+"\n")
		Rfich.write("mat=matrix(c(1,1,1,2),ncol=1)"+"\n")
		Rfich.write("layout(mat)"+"\n")
		Rfich.write("par(adj=0.5,bty=\"n\",font.lab=2,ps=10,las=1,lab=c(20,5,7),mar=c(2, 4, 4, 2)+0.1)"+"\n")
		# profil 1
		Rfich.write("barplot(mat1,border=NA,space=0,col=mycol,ylim=c(mn,mx),ylab=\"Reads per 10M reads\",xlab=\"\",main=\""+libName1+"\")"+"\n")
		Rfich.write("par(new=T)"+"\n")
		Rfich.write("barplot(mat2,border=NA,space=0,col=mycol,ylim=c(mn,mx))"+"\n")
		# pos
		Rfich.write("par(adj=0.5,bty=\"n\",font.lab=2,ps=10,las=1,lab=c(20,5,7),mar=c(5, 4, 2, 2)+0.1)"+"\n")
		Rfich.write("plot(seq("+st+","+sp+"-1,1),rep(0,"+str(int(sp)-int(st))+"),col=\"white\",type=\"h\",ylim=c(-1,1),yaxt=\"n\",main=\"\",xlab=\"Position (bp)\",ylab=\"\")"+"\n")
		Rfich.write("invisible(dev.off())"+"\n")
		Rfich.close()
		os.system("R --slave < "+myId+"_cmd.R")
	lign=fich.readline()

fich.close()
