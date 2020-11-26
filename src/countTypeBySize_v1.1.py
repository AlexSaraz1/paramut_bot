#!/usr/bin/python
# -*- coding: utf-8 -*-

import getopt,sys,os,string


usage='''

	Name :
	 		
	Synopsis :
	 		

'''




def readPriotity(allBankFile,bankName,bankPriority,allCount):
	allBankFich=open(allBankFile,'r')
	lign=allBankFich.readline()
	while (lign!=""):
		if (lign[0]!="#"):
			lignSplit=string.split(lign,sep="\t")
			bankName.append(lignSplit[0])
			bankPriority[lignSplit[0]]=int(lignSplit[1])
			allCount.append({})
		lign=allBankFich.readline()
	allBankFich.close()

def countTypeByLength(annotFile,bankPriority,allCount):
	fich=open(infile,'r')
	lign=fich.readline()
	while (lign != ""):
		if ((lign[0]!="#") and (lign[0:2]!="id")):
			lignSplit=string.split(lign,sep="\t")
			lg=int(lignSplit[2])
			count=int(lignSplit[3])
			typ=lignSplit[5]
			allCount[bankPriority[typ]-1][lg]+=count
		lign=fich.readline()
	fich.close()

if __name__=='__main__':
	## Lecture des info de BANQUES
	## Creation des listes de banques et dico de comptage 
	allBankFile=sys.argv[2]
	bankName=[]
	bankPriority={}
	allCount=[]
	readPriotity(allBankFile,bankName,bankPriority,allCount)
	## rajout d'un champ genomic ou non-annote
	bankName.append('nonAnnote')
	bankPriority['nonAnnote']=len(bankName)
	allCount.append({})
	## initialisation du dico de comptage a 0 pour toute les longueurs de séquence entre mini et maxi
	# longueur min/max
	mini=10
	maxi=51
	i=mini
	while (i < maxi):
		for aC in allCount:
			aC[i]=0
		i+=1	
	# lecture du fichier d'annotation et comptage des type en fonction de la longueur
	infile=sys.argv[1]
	countTypeByLength(infile,bankPriority,allCount)
	# output
	outName=sys.argv[3]
	outFich=open(outName,'w')
	i=mini
	# header
	outFich.write("lg")
	for b in bankName :
		outFich.write("\t"+b)
	outFich.write("\n")
	# count
	while ( i < maxi ):
		lg=i
		outFich.write(str(i))
		for bn in bankName :
			outFich.write("\t"+str(allCount[bankPriority[bn]-1][lg]))
		outFich.write("\n")
		i+=1
outFich.close()


