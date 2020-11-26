#!/usr/bin/python
# -*- coding: utf-8 -*-

import os,sys,string

pmodFile=sys.argv[1]
mappedTair10File=sys.argv[2]	
idLib=sys.argv[3]	


# reads .pmod file to get sequence id (idSeq) and reads count.
# informations are stored in the dictionary allSeq
# number of genomic position for each sequence is set at 0
allSeq={}
fich=open(pmodFile,'r')
lign=fich.readline()
while (lign != ""):
	idSeq=lign.replace("@","").replace("\n","")
	count=string.split(idSeq,sep="_")[1]
	lign=fich.readline()
	seq=lign.replace("\n","")
	lign=fich.readline()
	lign=fich.readline()
	if not allSeq.get(idSeq):
		allSeq[idSeq]=[seq,count,0]
	else:
		print "duplicated pmod "+idSeq
	lign=fich.readline()

fich.close()


# reads .sam file from bowtie2 and increment the number of genomic position for each sequence
fich=open(mappedTair10File,'r')
lign=fich.readline()
while (lign != ""):
	lignSplit=string.split(lign,sep="\t")
	if (lign[0] != "@"):
		if lignSplit[2]!="*":	
			IdSeq=lignSplit[0]		
			allSeq[IdSeq][2]+=1
	lign=fich.readline()

fich.close()



# reads .sam file from bowtie2 print gtf formated informations.
fich=open(mappedTair10File,'r')
lign=fich.readline()
while (lign != ""):
	lignSplit=string.split(lign,sep="\t")
	if (lign[0] != "@"):
		if lignSplit[2]!="*":	
			idSeq=lignSplit[0]		
			flg=lignSplit[1]
			st=lignSplit[3]
			seq=allSeq[idSeq][0]
			if flg=="0" or flg=="256":
				sd="+"
				print lignSplit[2]+"\tmySource\t"+idLib+"\t"+st+"\t"+str(int(st)+len(seq))+"\t.\t"+sd+"\t.\t"+idLib+" "+idLib+"_"+idSeq+" ; Count "+allSeq[idSeq][1]+" ; HitNumber "+str(allSeq[idSeq][2])+" ; Seq "+allSeq[idSeq][0]
			elif flg=="16" or flg=="272":
				sd="-"
				print lignSplit[2]+"\tmySource\t"+idLib+"\t"+st+"\t"+str(int(st)+len(seq))+"\t.\t"+sd+"\t.\t"+idLib+" "+idLib+"_"+idSeq+" ; Count "+allSeq[idSeq][1]+" ; HitNumber "+str(allSeq[idSeq][2])+" ; Seq "+allSeq[idSeq][0]
	lign=fich.readline()

fich.close()

