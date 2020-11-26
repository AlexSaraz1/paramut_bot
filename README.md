

# paramut_bot

---
##### Summary

This is a shell script that process a trimmed sRNA-seq fastq file.
It will:

1. Use processReads-median.pl from ncPRO-seq in order to regroup the reads with similar sequence
2. RUN bowtie2 in order to align the grouped reads on TAIR10 reference genome complemented with C-insert sequence as an additional chromosome. Output will be use for further processing
3. RUN bowtie2 in order to align the ungrouped reads on TAIR10 reference genome complemented with C-insert sequence as an additional chromosome. Output can be use for visualisation using IGV
4. Generate a gtf/gff file 
5. Annotate the sRNA sequences and generate an histogram of the proportion of reads versus reads length with annotation information
6. Generate a histogram comparing all sRNA reads length to the sRNA reads length of a given annotation (here C-insert)
7. Generate the distribution of sRNA reads over the C-insert sequence according to their length and strand.
8. Generate an analysis of the 5' nucleotide content of the sRNA reads aligned to the C-insert sequence.

---
##### Disclaimer

This pipeline has been made to help a colleage to process his sRNA sequencing data. It make use of several scripts that I wrote in python or R. This pipeline have been made for a specific type of data and processing and therefore might not be applicable for other purposes. 
This pipeline is not done to run on your personnal desktop.


##### Author

A. Sarazin

---
##### Dependancies

This pipeline has been tested on Ubuntu 15.04 with:

   * Python (v2.7.9)
   * Perl (v5.20.2)
   * R (3.6.3)
   * Bowtie2 (v2.4.1)
   * Samtools (v1.10)

---
##### Install/Usage

1. Download this directory and move inside
2. Edit the config.txt file 
2. Uncompress annotation files
```
gunzip Araport11_Annot/*.gff.gzip
```
3. Run
```
sh ./paramut_bot.sh <myFastqFile>
```

An example/test .fastq file, test1mReads.fastq, containing 1millions trimmed reads is included in the `./data` directory.

```
gunzip data/test1mReads.fastq.gz
sh ./paramut_bot.sh test1mReads.fastq
```

---
##### Configuration file

The file config.txt contains:
   * Path and name of the executable
   * Variables requiered for the processing
     * path to the bowtie2 index used for mapping        
     * path to a file containing annotation name, priority and gff formated files of annotations of interest (here annotations from Araport11 and one for Cinsert)
     * minimal and maximal length to consider for histograms.
     * Informations concerning the locus of interest for reads length histogram, sRNA profile and 5' nt content.
    

---
##### Workflow

A fastq file <myFastqFile> containing trimmed sRNA-seq data is place in `./data` directory.

The script run `processReads-median.pl` (taken from Chen *et al.* 2012 ; PMID:23044543) and generate a <myFastqFile>.pmod in `./data`. The .pmod file, is a fastq file were all reads with identical sequences are grouped and were the quality is the median quality of a given nucleotides in each of similar reads. Other files containing informations about read length and sequencing quality are also generated and store in `./data/groupReads`.

###### Mapping

Bowtie2 is call twice. Once on the newly created .pmod file and once on the initial .fastq file for visualisation using IGV. The bowtie2 index used for reads alignement is located in `./ref` directory. Here, TAIR10_chr_all_Cinsert bowtie2 index has been generated using TAIR10 reference after adding Cinsert sequence to the multifasta.
Mapped are redirected in `./mapping` directory.
For mapping of the ungrouped data, a bam file is generated and then sorted and indexed using samtools.
For mapping of grouped data, a .sam file is generated, used to create a gtf/gff file using `makeGff.py` and then deleted.
Note that both mapping are done with bowtie2 parameter `-k 500`, which is hard coded in the shell script.

The attribute field of the gtf/gff file contains informations about the number of reads of a given sequence, the number of genomic position as well as the sequence itself. The format of this file is important for the next processing steps.

###### Sequences annotation

This processed is done by comparing the genomic position of the sRNA-seq to various gff formated annotation files. Those annotations and path to the files are included in `../Araport11_Annot/Araport11Annot_Cinsert.txt`. This file also contain a priority value which is use for sequences with multiple genomic position or that overlap several annotation. If this is the case, the given sequence will be attributed the annotation with the lowest priority value (ie tRNA < rRNA < snoRNA ...). The annotations were taken from Araport11 and repeat-masker. The Cinsert.gff file contains a single line.

The script `./src/Annotate_v2.1.py` is used to attribute the annotation to the sRNA sequences. It require a fasta file in which the sequence header contains informations about the length, reads count and number of genomic position. This fasta file is generated by applying a awk command to the gff file and is store in the `./annotation` directory.
The output "annotated_"<libId>.txt is a table that contains one row per sequences with related informations, attributed annotation as well as identidier of the annotations it overlaps.

After this annotation procedure, a call to `./src/countTypeBySize_v1.1.py` is done in order to generate the table use to produce the annotation/length histograms. The table, call "countTypeLg_"<libId>.txt, is store in the `./annotation` directory and contain for each length, the number of reads for each annotation (or unannotated) for reads in the length range 10-51nt.
Table generated here serve as input for `./src/src/hist_size_annot.R` and `src/hist_size_annot_compare.R`.
    

###### sRNA profile

The sRNA profile over the C-insert sequence is generarated by `./src/profil4genes_1lib_stack_pos.py`. Argument for this script include `./profil/Cinsert.txt` containing the name of the locus of interrest, `./Araport11_Annot/Cinsert.gff` containing the informations of the locus of interrest, the gtf/gff file generated after the mapping procedure, the name of the sRNAlibrary and the number of mapped reads. The number of mapped reads is calculated using a awk command on the table generated by `./src/Annotate_v2.1.py`.

`./src/profil4genes_1lib_stack_pos.py`  will
   - apply a awk command to retrieve the sRNA nested over the locus of interrest
   - call `./src/sRNAdistri_plusStrand_gff3.py` to generated a table with the sRNA reads coverage of each nucleotide of the locus of interrest (one raw per nucleotide). The table contains 8 columns corresponding to the coverage for 20-21nt, 22-23nt, 24-25nt and other length for plus and minus strand.
   - write and call a R script to generated the graphical representation.

All files generated during this step are moved in `./profil/`.


###### sRNA 5' nucleotide content

This step make use of `./src/5primNucFromGff.py` on the gtf/gff file generated after mapping. It is call once without specification of genomic coordinate to generate a table with the 5' nucleotide content of all the reads and then once with specification of the genomic coordinate of the locus of interrest (Cinsert 1 5000) to generate the 5' nucleotide content of the reads from this locus. Any nucleotide different than A,T,G or C will be counted as N.

`./src/hist_5primNuc.R` is used to generate the graphical representation.
Files generated during this step are moved in `./5primNucleotide/`.


