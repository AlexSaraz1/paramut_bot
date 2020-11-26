#!/bin/sh


NORMAL="\\033[0;39m" 
RED="\\033[1;31m"
BLUE="\\033[0;34m"

. ./config.txt
echo "$RED""--- exec ---""$NORMAL"
# checking exec
listExec="$R_PATH $PERL_PATH $PYTHON_PATH $AWK_PATH $BOWTIE_PATH $SAMTOOLS_PATH"
for myEx in $listExec; do
if ! type "$myEx" > /dev/null; then
	echo "$RED"$myEx" not found""$NORMAL";
	echo ""
	exit 1
else
	echo $myEx;
fi
done
echo ""

echo ""
echo "$RED""--- options ---""$NORMAL"
echo "Bowtie2 index: "$BOWTIE2_INDEX
echo "Path to annotation files: "$ANNOT_BANK
echo "Small RNA length consider for annotation graph: "$MyMIN"-"$MyMAX
echo "Annotation use for reads length comparison: "$MyAnnot
echo "Annotation use for sRNA profil: "$MyAnnotForProfil
echo "GFF file use for sRNA profil: "$GffAnnotForProfil
echo "Genomic coordinate for 5 prime nucleotide comparison: "$MyGenomic_C":"$MyGenomic_St"-"$MyGenomic_Sp
echo ""
echo ""



# checking/creating directories
listVar="./mapping ./annotation ./profil ./5primNucleotide"
for myDir in $listVar; do
if [ ! -d "$myDir" ]; then
	echo $myDir" not found...creating...";
	mkdir $myDir;
fi
done
echo ""
echo ""

# loading module 
module load R
module load Bowtie2/2.0.2-goolf-1.4.10
module load SAMtools
echo "$RED""--- Running ---""$NORMAL"


FASTQ_FILE=$1
# Reads grouping
echo "-- Running processReads-median.pl (might take some times) --"
echo "$PERL_PATH ./src/processReads-median.pl -i $FASTQ_FILE -f phred33 -g 1 -D data -d data/groupReads/"
$PERL_PATH ./src/processReads-median.pl -i $FASTQ_FILE -f phred33 -g 1 -D data -d data/groupReads/
echo ""

echo "$RED""--- Mapping ---""$NORMAL"
## Mapping grouped reads
echo "-- Mapping grouped reads --"
#echo "$BOWTIE_PATH -p 14 -k 500 -x ref/TAIR10_chr_all_Cinsert -U data/$FASTQ_FILE.pmod --no-unal 2> mapping/$FASTQ_FILE.pmod.log > mapping/$FASTQ_FILE.pmod.sam"
#$BOWTIE_PATH -p 14 -k 500 -x ref/TAIR10_chr_all_Cinsert -U data/$FASTQ_FILE.pmod --no-unal 2> mapping/$FASTQ_FILE.pmod.log > mapping/$FASTQ_FILE.pmod.sam
echo "$BOWTIE_PATH -k 500 -x ref/TAIR10_chr_all_Cinsert -U data/$FASTQ_FILE.pmod --no-unal 2> mapping/$FASTQ_FILE.pmod.log > mapping/$FASTQ_FILE.pmod.sam"
$BOWTIE_PATH -k 500 -x ref/TAIR10_chr_all_Cinsert -U data/$FASTQ_FILE.pmod --no-unal 2> mapping/$FASTQ_FILE.pmod.log > mapping/$FASTQ_FILE.pmod.sam
echo ""

## Mapping ungrouped reads
echo "-- Mapping ungrouped reads for IGV visualisation --"
echo $FASTQ_FILE
#echo "$BOWTIE_PATH -p 14 -k 500 -x ref/TAIR10_chr_all_Cinsert -U data/$FASTQ_FILE --no-unal | samtools view -bSF4 - > mapping/$FASTQ_FILE.bam"
#$BOWTIE_PATH -p 14 -k 500 -x ref/TAIR10_chr_all_Cinsert -U data/$FASTQ_FILE --no-unal | $SAMTOOLS_PATH view -bSF4 - > mapping/$FASTQ_FILE.bam
echo "$BOWTIE_PATH -k 500 -x ref/TAIR10_chr_all_Cinsert -U data/$FASTQ_FILE --no-unal | samtools view -bSF4 - > mapping/$FASTQ_FILE.bam"
$BOWTIE_PATH -k 500 -x ref/TAIR10_chr_all_Cinsert -U data/$FASTQ_FILE --no-unal | $SAMTOOLS_PATH view -bSF4 - > mapping/$FASTQ_FILE.bam
echo ""

echo "-- Sorting mapped ungrouped reads and indexing --"
#
#echo "$SAMTOOLS_PATH sort mapping/$FASTQ_FILE.bam mapping/$FASTQ_FILE.sorted.bam"
#$SAMTOOLS_PATH sort mapping/$FASTQ_FILE.bam mapping/$FASTQ_FILE.sorted
#
echo "$SAMTOOLS_PATH sort mapping/$FASTQ_FILE.bam -o mapping/$FASTQ_FILE.sorted.bam"
$SAMTOOLS_PATH sort mapping/$FASTQ_FILE.bam -o mapping/$FASTQ_FILE.sorted.bam
echo "$SAMTOOLS_PATH index mapping/$FASTQ_FILE.sorted.bam"
$SAMTOOLS_PATH index mapping/$FASTQ_FILE.sorted.bam
rm mapping/$FASTQ_FILE.bam
echo ""


## create gff file for mapped grouped reads
libId=`  echo $FASTQ_FILE | sed 's/\.fastq//g' `
#echo $libId
gffFile=` echo "mapping/"$libId"_TAIR10.gff"`
#echo $gffFile
echo "-- Creating gff file for mapped grouped reads --"
echo "$PYTHON_PATH ./src/makeGff.py ./data/$FASTQ_FILE.pmod ./mapping/$FASTQ_FILE.pmod.sam $libId  > $gffFile"
$PYTHON_PATH ./src/makeGff.py ./data/$FASTQ_FILE.pmod ./mapping/$FASTQ_FILE.pmod.sam $libId  > $gffFile
rm ./mapping/$FASTQ_FILE.pmod.sam
echo ""


echo "$RED""--- Annotation ---""$NORMAL"
## Reads annotation
echo "-- Annotating mapped reads using $ANNOT_BANK --"
mappedFasta=` echo "annotation/"$libId"_seqQuiMatch.fa"`
AnnotData=` echo "annotation/annotated_"$libId".txt"`
$AWK_PATH '/^/{ print $10" "$12" "$13" "$15" "$16" "$NF}' $gffFile | sort -n | uniq | awk '/^/{ print ">"$1"\tlg="length($6)"\t"$2"="$3"\t"$4"="$5;print $6}' > $mappedFasta
echo "$PYTHON_PATH ./src/Annotate_v2.1.py -g $gffFile -f $mappedFasta -b $ANNOT_BANK -o $AnnotData"
$PYTHON_PATH ./src/Annotate_v2.1.py -g $gffFile -f $mappedFasta -b $ANNOT_BANK -o $AnnotData
echo ""

## Annotation graph
countLGFile=` echo "annotation/countTypeLg_"$libId".txt" `
#echo $countLGFile
echo "$PYTHON_PATH ./src/countTypeBySize_v1.1.py $AnnotData $ANNOT_BANK $countLGFile"
$PYTHON_PATH ./src/countTypeBySize_v1.1.py $AnnotData $ANNOT_BANK $countLGFile
$R_PATH --slave --args $ANNOT_BANK $countLGFile $libId $MyMIN $MyMAX < src/hist_size_annot.R
$R_PATH --slave --args $ANNOT_BANK $countLGFile $libId $MyMIN $MyMAX $MyAnnot < src/hist_size_annot_compare.R 
mv $libId".pdf" annotation/
mv $libId"_compare.pdf" annotation/
echo ""

echo "$RED""--- Profil and 5prime composition ---""$NORMAL"
## sRNA profil over Cinsert
echo "-- sRNA profil --"
NBreads=`$AWK_PATH 'BEGIN{ cmp=0} /^/{ if ($1 != "id"){cmp=cmp+$4}} END{print cmp}' $AnnotData ` 
### TODO : change calling to R and python in profil4genes_1lib_stack_pos.py
echo "$PYTHON_PATH ./src/profil4genes_1lib_stack_pos.py $MyAnnotForProfil $GffAnnotForProfil $gffFile $libId $NBreads"
$PYTHON_PATH ./src/profil4genes_1lib_stack_pos.py $MyAnnotForProfil $GffAnnotForProfil $gffFile $libId $NBreads
mv distriSur*.pdf profil/
mv ReadLength* profil/
mv reads*.gff profil/
mv nbReadsStranded* profil/
mv *_cmd.R profil/
mv *_cmd.awk profil/
echo ""

## 5 prime nucleotide composition
echo "-- 5 prime nucleotide composition --"
FiveprimNameTotal=` echo $libId"_allReads" `
echo "$PYTHON_PATH ./src/5primNucFromGff.py -g $gffFile -o $FiveprimNameTotal"
$PYTHON_PATH ./src/5primNucFromGff.py -g $gffFile -o $FiveprimNameTotal 
$R_PATH  --slave --args $FiveprimNameTotal"_read_length_5primNuc.txt" $FiveprimNameTotal $MyMIN $MyMAX < ./src/hist_5primNuc.R 
echo "--"
FiveprimNameLocus=` echo $libId"_"$MyGenomic_C"_"$MyGenomic_St"_"$MyGenomic_Sp`
echo "$PYTHON_PATH ./src/5primNucFromGff.py -g $gffFile -o $FiveprimNameLocus -c $MyGenomic_C -b $MyGenomic_St -e $MyGenomic_Sp"
$PYTHON_PATH ./src/5primNucFromGff.py -g $gffFile -o $FiveprimNameLocus -c $MyGenomic_C -b $MyGenomic_St -e $MyGenomic_Sp
$R_PATH  --slave --args $FiveprimNameLocus"_read_length_5primNuc.txt" $FiveprimNameLocus $MyMIN $MyMAX < ./src/hist_5primNuc.R 
mv "proportion_test1mRead"*.pdf 5primNucleotide/
mv $libId"_allReads.pdf" 5primNucleotide/
mv $FiveprimNameLocus*".pdf" 5primNucleotide/
mv $FiveprimNameLocus*".txt" 5primNucleotide/
mv *read_length_5primNuc* 5primNucleotide/
echo ""

