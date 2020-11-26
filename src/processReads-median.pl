###!/usr/bin/perl -w
#created by Chongjian CHEN- 16/09/2011

#The script is used to check the input format, calcualte read length distribution, compute median quality score of each base
#if -g is specified with 1, group reads which have the same sequence in raw sequencing data. For fastq reads, the median quality score is computed in each read group.
#or else, just output original sequence with extended read name with "_1" suffix in order to have the same name format as group read.

use Getopt::Std;
use utf8;
use strict;

## get options from command line
my %opts= ('i'=>'','f'=>'','g'=>'','D'=>'','d'=>'');

sub usage{
print STDERR <<EOF;
    usage:  $0 -i read_file -f format -g 1|0 -D seq_dir -d data_dir [-h]

     -h   : help message;
     -i   : input read file;
     -f   : format of read file, could be specified as "fa","csfast", "solexa", "solexa1.3", and "phred33";
     -g   : group read or not, 1: yes, 0: no;
     -D   : the raw sequence directory;
     -d   : the data directory;

    example1: $0 -i file -f "phred33" -g 1 -D "rawdata" -d "data"
    example2: $0 -i file -f "fa" -g 1 -D "rawdata" -d "data"
    example3: $0 -i file -f "csfast" -g 1 -D "rawdata" -d "data"

EOF
    exit;
}

getopts('i:f:g:D:d:h', \%opts) || usage();
usage() if ($opts{h});

my ($read_file,$format,$do_group,$seq_dir,$data_dir)=($opts{i},$opts{f},$opts{g},$opts{D},$opts{d});

##specify start identifier of read (> or @) in sequence file
my $stid=">";
if($format!~/fa/){
    $stid="@";
}

my ($max_len,$cur_seq,$cur_qc)=(0,"","");
my %seq_len=(); # hash of sequence length
my %seq_group=(); # hash of read group 
my %seq_qc=(); # hash of quality score
my %seq_base=(); # hash of read base 

open (SEQ,"$seq_dir/$read_file") || warn "can't open $read_file";
open (GROUPOUT,">$seq_dir/$read_file.pmod"); #modified read file after processing
while(<SEQ>){
    if(/^$stid(\S+)/){
	my $cur_header=$1;
	$cur_seq=<SEQ>;
	chomp($cur_seq);
	my $nuc_pos=1;
	foreach my $c (split(//,$cur_seq)){
	    $seq_base{$nuc_pos}{$c}++;
	    $nuc_pos++;
}
	my $cur_seq_len=length($cur_seq);
        ##validate the fasta/csfast format
	if($format eq "fasta"){
	    if($cur_seq!~/[ATGCNatgcN]{$cur_seq_len,}/){
		`rm -f $seq_dir/$read_file.group`;
		die("$seq_dir/$read_file is not a fasta file\n");
}
            #output read without grouping
	    if(!$do_group){
		print GROUPOUT ">",$cur_header,"_1\n",$cur_seq,"\n";
}
}
	if($format eq "csfasta"){
	    my $len_nonuc=$cur_seq_len-1;
	    ##Remove the first tag from the sequence length
	    $cur_seq_len=$cur_seq_len-1;
	    if($cur_seq!~/[ATGCatgc][0123\.]{$len_nonuc,}/){
		`rm -f $seq_dir/$read_file.group`;
		die("$seq_dir/$read_file is not a csfasta file\n");
}
            #output read without grouping
	    if(!$do_group){
		print GROUPOUT ">",$cur_header,"_1\n",$cur_seq,"\n";
}
}
	#get length info
	$max_len=$cur_seq_len if ($cur_seq_len>$max_len);
	$seq_len{$cur_seq_len}++;
	#get read group (distinct read) information
	$seq_group{$cur_seq}{"ct"}++;
	#process fastq reads
	if($stid eq "@"){
	    my $mid_line=<SEQ>;
            #validate the fastq format
	    if($mid_line!~/^\+/){
		`rm -f $seq_dir/$read_file.group`;
		die("$seq_dir/$read_file is not a fastq file\n");
}
	    $cur_qc=<SEQ>;
	    #output read without grouping
	    if(!$do_group){
		print GROUPOUT "@",$cur_header,"_1\n",$cur_seq,"\n","+",$cur_header,"_1\n",$cur_qc;
}
	    chomp($cur_qc);
	    my $pos=1;
	    foreach my $c (split(//, $cur_qc)){
		#get unicode value from characters
		my $cur_ord=ord($c);
		#validate the fastq format
		if(($format eq "phred33") && ($cur_ord<33 || $cur_ord>80)){
		    `rm -f $seq_dir/$read_file.group`;
		    die("$seq_dir/$read_file is not a sanger fastq file\n");
}
		if(($format eq "solexa") && ($cur_ord<59 || $cur_ord>126)){
		    `rm -f $seq_dir/$read_file.group`;
		    die("$seq_dir/$read_file is not a solexa 1.0 fastq file\n");
}
		if(($format eq "solexa1.3") && ($cur_ord<64 || $cur_ord>104)){
#		    warn $cur_ord,"\n";
		    `rm -f $seq_dir/$read_file.group`;
		    die("$seq_dir/$read_file is not a solexa 1.3 fastq file\n");
}
		#store unicode value of positional quality score
		if (!$seq_qc{$pos}){
		    $seq_qc{$pos}=$cur_ord;
}
		else{
		    $seq_qc{$pos}.="\t" . $cur_ord;
}
		#store base quality score info for each read group, if $do_group is active
		if($do_group==1){
		    if(!$seq_group{$cur_seq}{$pos}){
			$seq_group{$cur_seq}{$pos}=$cur_ord;
}
		    else{
			$seq_group{$cur_seq}{$pos}.="\t" . $cur_ord;
}
}
		$pos++;
}
}
}
}
close(SEQ);

##get sample name
my $sample="";
if($read_file=~/(\S+)\./){
    $sample=$1;
}

##output length distribution of abundant reads
open (LENOUT,">$data_dir/$sample\_readlen.data");
print LENOUT "idx\t$sample\n";
foreach my $len (sort {$a <=> $b} keys %seq_len){
    print LENOUT $len,"\t",$seq_len{$len},"\n";
}
close(LENOUT);

##output read base composition 
open (BASE,">$data_dir/$sample\_basestat.data");
open (GC,">$data_dir/$sample\_baseGCstat.data");
print BASE "idx\tA\tT\tG\tC\n";
print GC "idx\t$sample\n";
foreach my $pos (1..$max_len){
    $seq_base{$pos}{"A"}=0 if (!$seq_base{$pos}{"A"});
    $seq_base{$pos}{"T"}=0 if (!$seq_base{$pos}{"T"});
    $seq_base{$pos}{"G"}=0 if (!$seq_base{$pos}{"G"});
    $seq_base{$pos}{"C"}=0 if (!$seq_base{$pos}{"C"});
    my $cur_basect=$seq_base{$pos}{"A"}+$seq_base{$pos}{"T"}+$seq_base{$pos}{"G"}+$seq_base{$pos}{"C"};
    my $cur_gcct=$seq_base{$pos}{"G"}+$seq_base{$pos}{"C"};
    print BASE $pos,"\t",100*$seq_base{$pos}{"A"}/$cur_basect,"\t",100*$seq_base{$pos}{"T"}/$cur_basect,"\t",100*$seq_base{$pos}{"G"}/$cur_basect,"\t",100*$seq_base{$pos}{"C"}/$cur_basect,"\n";
    print GC $pos,"\t",100*$cur_gcct/$cur_basect,"\n";
}
close(BASE);
close(GC);

##output length distribution of distinct reads
open (DISTLENOUT,">$data_dir/$sample\_distinct_readlen.data");
print DISTLENOUT "idx\t$sample\n";
my %distinct_len;
foreach my $cur_read (keys %seq_group){
    my $cur_len=length($cur_read);
    $distinct_len{$cur_len}++;
}
foreach my $len  (sort {$a <=> $b} keys %distinct_len){
    print DISTLENOUT $len,"\t",$distinct_len{$len},"\n";
}
close(DISTLENOUT);

##output read groups  if $do_group is active
if($do_group==1){
    my $gp_id=1;
    foreach my $sg (keys %seq_group){
        #for fastq format read
	if($stid eq "@"){
	    #use group id and the number of read in the group to uniquely define a read group
	    print GROUPOUT "\@$gp_id\_",$seq_group{$sg}{"ct"},"\n";
	    print GROUPOUT $sg,"\n";
	    print GROUPOUT "\+$gp_id\_",$seq_group{$sg}{"ct"},"\n";
	    my $cur_seq_len=length($sg);
	    foreach my $id (1..$cur_seq_len){
		my $cur_median_qc=median_qc($seq_group{$sg}{$id});
		print GROUPOUT chr($cur_median_qc);
}
	    print GROUPOUT "\n";
	    $gp_id++;
}
        #for fa and csfast reads
	else{ 
	    #use group id and the number of read in the group to uniquely define a read group
	    print GROUPOUT ">$gp_id\_",$seq_group{$sg}{"ct"},"\n";
	    print GROUPOUT $sg,"\n";
	    $gp_id++;
}
}
}
close(GROUPOUT);

#output median quality score for each position in all fastq reads
if($stid eq "@"){
    open(QCOUT,">$data_dir/$sample\_medquality.data");
    print QCOUT "idx\t$sample\n";
    foreach my $pos (1..$max_len){
	my $cur_qcstr=$seq_qc{$pos};
	my $pos_mediansc=median_qc($cur_qcstr);
	$pos_mediansc=qc_real_sc($pos_mediansc);
	print QCOUT $pos,"\t",$pos_mediansc,"\n";
}
    close(QCOUT);
}

#########
#function to get real quanlity score based on the different fastq format
sub qc_real_sc{
    my ($sc)=@_;
    if($format eq "solexa1.3"){
	$sc-=64;
}
    elsif($format eq "phred33"){
	$sc-=33;
}
    else{
	$sc-=64;
	$sc=log(1+10^($sc/10))/log(10);
}
}

#########
#function to calcluate the median score in a given string.
sub median_qc{
    my ($s)=@_;
    my @array=split(/\t/,$s);
    my $count=$#array + 1;
    @array=sort {$a <=> $b} @array;
    my $median="";
    if ($count % 2){
	$median=$array[int($count/2)];
} 
    else{
	$median=($array[$count/2]+$array[$count/2 - 1])/2;
}
    return($median);
}
