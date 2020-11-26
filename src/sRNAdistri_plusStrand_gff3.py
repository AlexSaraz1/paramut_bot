#!/usr/bin/python


import sys,os,string

TEname=sys.argv[2]
TEst=int(sys.argv[3])
TEsp=int(sys.argv[4])
TElg=TEsp-TEst
sc21plus=[]
sc22plus=[]
sc24plus=[]
sc21minus=[]
sc22minus=[]
sc24minus=[]
scOtherplus=[]
scOtherminus=[]
i=1
while i < TElg+1 :
	sc21plus.append(0)
	sc22plus.append(0)
	sc24plus.append(0)
	sc21minus.append(0)
	sc22minus.append(0)
	sc24minus.append(0)
	scOtherplus.append(0)
	scOtherminus.append(0)
	i+=1

infile=sys.argv[1]
fich=open(infile,'r')
lign=fich.readline()
while lign != "" :
	lignSplit=string.split(lign.replace("\n",""),sep="\t")
	sChro=lignSplit[0]
	sSt=int(lignSplit[3])
	sSp=int(lignSplit[4])
	sSd=lignSplit[6]
	attribute=lignSplit[8]
	#print attribute
	attributeSplit=string.split(attribute.replace(" ; "," "),sep=" ")
	sCount=float(attributeSplit[3])
	sSeq=attributeSplit[7]
	sLg=len(sSeq)
	sSt_C=sSt-TEst
	sSp_C=sSp-TEst
	if (sSd=="+"):
		if (sLg==21 or sLg==20) :
			i=0
			while i < sLg :
				sc21plus[sSt_C+i]=sc21plus[sSt_C+i]+sCount
				i+=1
		elif (sLg==22 or sLg==23) :
			i=0
			while i < sLg :
				sc22plus[sSt_C+i]=sc22plus[sSt_C+i]+sCount
				i+=1
		elif (sLg==24 or sLg==25) :
			i=0
			while i < sLg :
				sc24plus[sSt_C+i]=sc24plus[sSt_C+i]+sCount
				i+=1
		else :
			i=0
			while i < sLg :
				scOtherplus[sSt_C+i]=scOtherplus[sSt_C+i]+sCount
				i+=1
	elif (sSd=="-") :
		if (sLg==21 or sLg==20) :
			i=0
			while i < sLg :
				sc21minus[sSt_C+i]=sc21minus[sSt_C+i]+sCount
				i+=1
		elif (sLg==22 or sLg==23) :
			i=0
			while i < sLg :
				sc22minus[sSt_C+i]=sc22minus[sSt_C+i]+sCount
				i+=1
		elif (sLg==24 or sLg==25) :
			i=0
			while i < sLg :
				sc24minus[sSt_C+i]=sc24minus[sSt_C+i]+sCount
				i+=1
		else :
			i=0
			while i < sLg :
				scOtherminus[sSt_C+i]=scOtherminus[sSt_C+i]+sCount
				i+=1
	lign=fich.readline()

fich.close()

n=0
while n < len(sc21plus) :
	print str(sc21plus[n])+"\t"+str(sc22plus[n])+"\t"+str(sc24plus[n])+"\t"+str(scOtherplus[n])+"\t"+str(sc21minus[n])+"\t"+str(sc22minus[n])+"\t"+str(sc24minus[n])+"\t"+str(scOtherminus[n])
	n+=1

