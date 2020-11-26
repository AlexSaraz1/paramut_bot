#!/usr/bin/python


import getopt,sys,os,string

usage='''

	Name :
	 		Annotate_v2.0.py
	Synopsis :
	 		./Annotate_v2.0.py -g <gff_file> -f <fasta_file> -b <bank_file> -o outfile
	
	Description :
	 	Annotate compare les positions des features present dans le fichier .gff <gff_file> avec les annotations des fichiers .gff reference dans le fichier <bank_file>.
		Un feature sera annote si il est entierement inclus dans l'une des annotation de reference.
	
	Options :
	 	-g, --gffFile
			feature file in gff format (REQUIRE)
		
		-f, --fastaFile
			feature file in fasta format (REQUIRE)
		
		-b, --banksFile
			liste of reference annotations file (REQUIRE)
			format : 
	   			Bank_One_Name	number_one	bank_one_file_pathway
				Bank_two_Name	number_two	bank_two_file_pathway
				(...)
		-o, --outfile
			output filename

	A.S. 13/04/2010

'''

## lecture fichier gff et stockage dans une matrice	
def readGff(namefile,mset):
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
			if (mset.get(chro)):
				mset[chro].append([chro,src,typ,st,sp,scr,sens,phase,ident])
			else :
				mset[chro]=[]
				mset[chro].append([chro,src,typ,st,sp,scr,sens,phase,ident])
		lign=fichSet.readline()
	print "read file "+namefile+" Completed"
	fichSet.close()
	return mset



## fonction d'annotation 
# mset1 = liste de feature a annoter (matchs)
# mset2 = liste d'annotation
# dicoLst pour le stockage des annotations (regions)
# bankPriority = priorite des annotations (regions)
def annote(mset1,mset2,dicoLst,bankPriority):
	# limite des compteurs pour affichage de l'avance de l'annotation
	limA=100 # limite du compteur d'avance dans le region ; qd atteinte on affiche
	limM=10000
	# compteurs pour l'annotation
	m=0 # compteur du nombre de matchs
	m2=0 # compteur du nombre de region
	flg=0 # flag sur le nombre de match
	flg2=0
	# pour chaque region
	for r in mset2 :
		# on reprend la ou etait le flag du nombre de match
		m=flg
		# incrementation du compteur de region
		m2+=1
		# recuperation des info de region
		chroR=r[0]
		typ=r[2]
		prioNb=bankPriority[typ]-1
		stR=int(r[3])
		spR=int(r[4])
		idR=r[8]
		# tant qu'on est pas au bout de la liste de matchs 
		while ( m < len(mset1) ):
			# recuperation des info de region
			chroM=mset1[m][0]
			stM=int(mset1[m][3])
			spM=int(mset1[m][4])
			#sc=int(mset1[m][5])
			sc=mset1[m][5]
			lg=int(spM)-int(stM)+1
			idM=mset1[m][8]
			# sRNA avant region
			if (spM <= stR ) :
				m+=1
			# sRNA overlap le debut de la region
			elif ( (stM < stR) and (spM >= stR) ):
				# c'est le 1er a overlapper, on le flag
				if (flg2==0):
					flg=m
					flg2=1
				# c'est pas le 1er, on le flag pas
				else : m+=1
			# match inclu dans region
			elif (( stM >= stR ) and (spM <= spR)):
				# "good one" on stock 
				# verification des id de region attribue a ce match pour eviter la redondance
				# si ce match a deja ete annote
				if (dicoLst[prioNb].get(idM)):
					ok=0 # flag de redondance
					# on verifie les id de regions pour eviter la redondance
					for t in dicoLst[prioNb][idM]:
						# si id de region deja presente
						if(t==idR):
							ok=1
							pass
					# si flag de redondance = 0 on est pas redondant et on ajoute 
					if (ok==0):
						dicoLst[prioNb][idM].append(idR)
				# si aucune annotation pour ce match pas de risque de redondance
				else :
					dicoLst[prioNb][idM]=[]
					dicoLst[prioNb][idM].append(idR)
				# on incremente le compteur de match
				m+=1
			# si on n'a pas fini les regions mais que l'on a fini les matchs
			## facultatif?!
			elif (m2<len(mset2) and m==(len(mset1)-1)) :
				m=flg
				flg2=0
				break
			# si match apres region
			elif (stM > spR):
				m=flg
				flg2=0
				break
			# 
			else :
				m+=1
	return dicoLst


if __name__=='__main__':
	# recup des options
	try:
		opts,args = getopt.getopt(sys.argv[1:], "ho:g:f:b:o", ["help", "gffFile=","fastaFile=","banksFile=","outFile="])
	except getopt.GetoptError, err:
		# print help information and exit:
		print str(err) # will print something like "option -a not recognized"
		sys.stderr.write(usage)
		sys.exit(2)
	# fichier des match a annoter .gff 
	INFILE=""
	# fichier des match a annoter .fa
	fastaFile=""
	# fichier contenant les infi de bank (nom, priorite, fichier)
	banksFile=""
	outFile=""
	for o, a in opts:
		if o == "-v":
			verbose = True
		elif o in ("-h","--help"):
			print usage
			sys.exit()
		elif o in ("-g","--gffFile"):
			INFILE = a
		elif o in ("-f", "--fastaFile"):
			fastaFile = a
		elif o in ("-b", "--banksFile"):
			banksFile = a
		elif o in ("-o", "--outFile"):
			outFile = a
		else:
			assert False, "unhandled option"
	print "#"+INFILE+" "+fastaFile+" "+banksFile
	# si pas de fichier en input/output
	if ((INFILE=="") or (fastaFile=="") or (banksFile=="") or (outFile=="")):
		sys.stderr.write(usage)
		sys.exit()
	# lecture fichier 1 (gff) et stockage dans un dico pour chque chromosome	
	mset1={}
	set1Name=INFILE
	mset1=readGff(set1Name,mset1)
	# trie pour chaque chromosome en fonction du start et stop des regions
	for c in mset1.keys() :
		mset1[c].sort(key=lambda x:x[4])
		mset1[c].sort(key=lambda x:x[3])
	## BANQUES
	bankFile=[]
	bankName=[]
	bankPriority={}
	allAnnote=[]
	# lecture du fichier contenant les info des banks
	allBankFile=banksFile
	allBankFich=open(allBankFile,'r')
	lign=allBankFich.readline()
	while (lign!=""):
		if (lign[0]!="#"):
			lignSplit=string.split(lign,sep="\t")
			bankName.append(lignSplit[0])
			if bankPriority.get(lignSplit[0]):
				print "duplicated type "+lignSplit[0]+" in bankFile "+banksFile
			else :
				bankPriority[lignSplit[0]]=int(lignSplit[1])
			bankFile.append(lignSplit[2].replace("\n",""))
			allAnnote.append({})
		lign=allBankFich.readline()
	allBankFich.close()
	## lecture de tout les fichiers d'annotation et stockage dans une matrice
	i=0
	mset2={}
	while (i<len(bankName)):
		## lecture fichier gff et stockage dans une matrice
		set2Name=bankFile[i]
		mset2=readGff(set2Name,mset2)
		i+=1
	# trie pour chaque chromosome en fonction du start et stop des regions
	for c in mset2.keys() :
		mset2[c].sort(key=lambda x:x[4])
		mset2[c].sort(key=lambda x:x[3])
	# annote chromosome par chromosome
	for c in mset2.keys():
		if mset1.get(c):
			print "----"
			print "chromosome "+c
			annote(mset1[c],mset2[c],allAnnote,bankPriority)
	##
	fastaName=fastaFile
	fastaFich=open(fastaName,'r')
	names={}
	hits={}
	counts={}
	seqs={}
	print "----"
	print "lecture du fasta"
	lign=fastaFich.readline()
	# compteur de seqeunces
	nbSeq=0
	nbSeqUniq=0
	# pour toute les lignes du fasta
	while (lign != "") :
		# si c'est le header de la sequence
		if ( lign[0]==">" ) :
			# incrementation du compteur de sequences
			lignsplit=string.split(lign.replace("\n","").replace(">",""),sep="\t")
			idseq=lignsplit[0]
			# verification que l'on a pas 2 fois le mm identifiant
			if ( names.get(idseq) ):
				print "2 sequences have the same id :"+seqName
				print "only one will be use."
			# sinon retient l'id
 			else :
				names[idseq]=idseq
				counts[idseq]=lignsplit[2].replace("Count=","")
				hits[idseq]=lignsplit[3].replace("HitNumber=","")
				nbSeqUniq+=1
				nbSeq=nbSeq+int(lignsplit[2].replace("Count=",""))
		# si c'est pas le header c'est la sequence
		else :
			# on choppe la sequence
			if ( seqs.get(idseq) ):
				seqlst=[seqs[idseq],lign.replace("\n","")]
				seqs[idseq]=string.join(seqlst,sep="")	
			else :
				seqs[idseq]=[lign.replace("\n","")]
		lign=fastaFich.readline()
	
	fastaFich.close()
	print "fasta finish"
	## tableau regroupant les annotations
	outfich=open(outFile,'w')
	outfich.write("id\tsequence\tlength\tnumber of de seq\tnumber of hit\tAnnotation\t"+string.join(bankName,sep="\t")+"\n")
	#print "id\tsequence\tlength\tnumber of de seq\tnumber of hit\tAnnotation\t"+string.join(bankName,sep="\t")
	for nam in names.keys():
		lignOut=[]
		lignOut.append(nam)
		chmpSeq=seqs[nam][0]
		lignOut.append(chmpSeq)
		lignOut.append(str(len(chmpSeq)))
		lignOut.append(counts[nam])
		lignOut.append(hits[nam])
		annote=""
		lignOutbis=[]
		i=0
		while (i < len(bankName)):
			j=bankPriority[bankName[i]]-1
			if allAnnote[j].get(nam):
				lignOutbis.append(string.join(allAnnote[j][nam]," ; "))
				if (annote==""):
					annote=bankName[i]
			else :
				lignOutbis.append("-")
			i+=1
		if (annote==""):
			annote="nonAnnote"
		outfich.write(string.join(lignOut,sep="\t")+"\t"+annote+"\t"+string.join(lignOutbis,sep="\t")+"\n")
	outfich.close()
	del(mset1)
	del(mset2)
	del(allAnnote)


