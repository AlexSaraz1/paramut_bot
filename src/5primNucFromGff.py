#!/usr/bin/python
# -*- coding: utf-8 -*-


import getopt,sys,os,string

usage='''

	Name :
	 		5primNucFromGff.py
	Synopsis :
	
	Description :
	
	Options :
	 	-g, --gffFile
			feature file in gff format (REQUIRE)
		
		-o, --outname
			prefix of output files (REQUIRE)

		-c, --chromosome
			chromosome of the selected region
		-b, --begin
			beginning of the selected region
		-e, --end
			end of the selected region

	A.S. 08/2017

'''
## lecture fichier gff et stockage des sequences, id et reads count dans un dico
def readGff(namefile,allSeq):
	#mset=[]
	setName=namefile
	fichSet=open(setName,'r')
	lign=fichSet.readline()
	print "----"
	print "reading file "+namefile
	while ( lign != "" ):
		if (lign[0]!="#"):
			lignSplit=string.split((lign.replace("="," ")).replace("\n",""),"\t")
			chro=lignSplit[0]
			src=lignSplit[1]
			typ=lignSplit[2]
			st=int(lignSplit[3])
			sp=int(lignSplit[4])
			scr=lignSplit[5]
			sens=lignSplit[6]
			phase=lignSplit[7]
			group=lignSplit[8]
			# split et lecture de l'id du feature
			ident=string.split(group.replace("="," ").replace(";"," ; ").replace("\"","")," ")[1]
			seq=string.split(group.replace("="," ").replace(";"," ; ").replace("\"","")," ")[16]
			nbR=string.split(group.replace("="," ").replace(";"," ; ").replace("\"","")," ")[6]
			if not allSeq.get(seq):
				allSeq[seq]=[seq,ident,nbR]
		lign=fichSet.readline()
	print "read file "+namefile+" Completed"
	fichSet.close()
	return allSeq

## lecture fichier gff et stockage des sequences, id et reads count dans un dico
## filtre en fonction de chromosome begin end
def readGffFilter(namefile,allSeq,chromosome,begin,end):
	#mset=[]
	setName=namefile
	fichSet=open(setName,'r')
	lign=fichSet.readline()
	print "--"
	print "reading file "+namefile
	while ( lign != "" ):
		if (lign[0]!="#"):
			lignSplit=string.split((lign.replace("="," ")).replace("\n",""),"\t")
			chro=lignSplit[0]
			src=lignSplit[1]
			typ=lignSplit[2]
			st=int(lignSplit[3])
			sp=int(lignSplit[4])
			scr=lignSplit[5]
			sens=lignSplit[6]
			phase=lignSplit[7]
			group=lignSplit[8]
			# split et lecture de l'id du feature
			ident=string.split(group.replace("="," ").replace(";"," ; ").replace("\"","")," ")[1]
			seq=string.split(group.replace("="," ").replace(";"," ; ").replace("\"","")," ")[16]
			nbR=string.split(group.replace("="," ").replace(";"," ; ").replace("\"","")," ")[6]
			if chro==chromosome and st >= begin and sp <= end:
				if not allSeq.get(seq):
					allSeq[seq]=[seq,ident,nbR]
		lign=fichSet.readline()
	print "read file "+namefile+" Completed"
	fichSet.close()
	return allSeq


# dico that will contain count, length, 5 prime Nt
def createDico(dico,lgMin,lgMax):
	i=lgMin
	while i <= lgMax :
		dico[i]=[0,0,0,0,0]
		i+=1
	return dico


if __name__=='__main__':
	# recup des options
	try:
		opts,args = getopt.getopt(sys.argv[1:], "ho:g:c:b:e:o", ["help", "gffFile=","outname=","chromosome=","begin=","end="])
	except getopt.GetoptError, err:
		# print help information and exit:
		print str(err) 
		sys.stderr.write(usage)
		sys.exit(2)
	# fichier des match a annoter .gff 
	INFILE=""
	# prefix of output files (REQUIRE)
	outname=""
	#
	chromosome=""
	begin=""
	end=""
	for o, a in opts:
		if o == "-v":
			verbose = True
		elif o in ("-h","--help"):
			print usage
			sys.exit()
		elif o in ("-g","--gffFile"):
			INFILE = a	
		elif o in ("-o","--outname"):
			outname = a	
		elif o in ("-c", "--chromosome"):
			chromosome = a
		elif o in ("-b","--begin"):
			begin = int(a)
		elif o in ("-e", "--end"):
			end = int(a)
		else:
			assert False, "unhandled option"
	# si pas de fichier en input/output
	if ((INFILE=="") or (outname=="")):
		sys.stderr.write(usage)
		sys.exit()
	
	# dico 
	lgMin=18
	lgMax=26
	allNuc={}
	allNuc=createDico(allNuc,lgMin,lgMax)
	
	#
	allSeq={}
	if chromosome=="" or begin=="" or end == "" :
		print "--"
		print "chromosome/begin/end not define"
		print "counting 5prime nucleotide composition for all reads"
		allSeq=readGff(INFILE,allSeq)
	else :
		allSeq=readGffFilter(INFILE,allSeq,chromosome,begin,end)
	
	#
	nucDico={"A":0,"T":1,"G":2,"C":3,"N":4}
	# count 5 prim nt and fill allNuc
	for k in allSeq.keys():
		myNuc=k[0]
		nbR=int(allSeq[k][2])
		lg=len(k)
		nucNb=nucDico[myNuc]
		allNuc[lg][nucNb]=allNuc[lg][nucNb]+nbR
	
	# return dico
	outfich=open(outname+"_read_length_5primNuc.txt",'w')
	i=lgMin
	outfich.write(string.join(nucDico.keys(),sep="\t")+"\n")
	while i <= lgMax:
		tmp=[]
		for n in nucDico.keys():
			mycount=allNuc[i][nucDico[n]]
			tmp.append(str(mycount))
		outfich.write(string.join(tmp,sep="\t")+"\n")		
		i+=1
	#
	outfich.close()

