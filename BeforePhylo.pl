#!/usr/bin/perl -w

# BeforePhylo (version 0.9.0): a Perl script for manipulating multiple sequence alignments for phylogenetic reconstruction
# Copyright (C) 2012-2014, Qiyun Zhu. All rights reserved.
# Licensed under BSD 2-clause license.

use warnings;
use strict;
$| = 1;

## Print program information ##

print "
Usage:
  perl BeforePhylo.pl [option(s)] [MSA(s)]

Input:
  One or more trees multiple sequence alignment (MSA) files in FASTA format.

Options:
  -type=(dna|protein|condon): type of data (default=dna)
  -filter=(taxa list file): only retain sequences defined in a taxa list file
  -conc=(none|raxml|beast|mrbayes): concantenate multiple MSAs and generate partition table
  -sort: sort sequence names in alphabetical order
  -Gblocks=(Gblocks program): perform Gblocks on each MSA
  -trim: remove empty sites
  -fillends: replace initial and final gaps with 'N's
  -fillgaps=(cutoff, default=10): replace in-sequence gaps no shorter than the cutoff with 'N's
  -fill: fill both ends and gaps
  -unalign: remove all gaps (make the sequences unaligned)
  -N: replace ambiguous codes with 'N's
  -10: keep the first 10 characters of sequence names
  -codon: divide dataset into three codon positions
  -numerize: translate taxon names into numbers, to avoid too complicated / duplicated names
  -translate=(dictionary file): translate taxon names according to a dictionary
  -partition=(RAxML-style partition file): split a master alignment by partition
  -output=(fasta|nexus|phylip|pir): output file format (default=fasta)
  -overwrite: overwrite original files (default is to create new files)

Examples:
  perl BeforePhylo.pl -type=codon -trim -fill -Gblocks coI.fas coII.fas
  perl BeforePhylo.pl -filter=myTaxa.txt -output=nexus -unalign *.fas
  perl BeforePhylo.pl -conc=mrbayes *.fas
\n"
and exit 0 if ($#ARGV < 0);


## Global variables ##

my $s; my $t;				# buffer string
my $i; my $j;				# buffer integer
my @a;						# buffer array
my %h;						# buffer hash

my $key;					# buffer key

my @files;					# partition file names
my %seqs;					# names of the taxa with sequences
my %parts;					# partition names vs sizes
my @part_order;				# partition names in order


## Define arguments ##

my %args =(
	"type" => "dna",
	"format" => "fasta",
	"filter" => 0,
	"concantenation" => 0,
	"Gblocks" => 0,
	"overwrite" => 0,
	"sort" => 0,
	"enize" => 0,
	"fillends" => 0,
	"fillgaps" => 0,
	"deEmpty" => 0,
	"unalign" => 0,
	"codon" => 0,
	"numerize" => 0,
	"translate" => 0,
	"partition" => 0
);

## read arguments ##

for ($i=0; $i<=$#ARGV; $i++){
	for ($ARGV[$i]){
		if (/^-/){
			/^-type=(.+)$/ and $args{"type"} = $1;
			/^-filter=(.+)$/ and $args{"filter"} = $1;
			/^-conc/ and $args{"concantenation"} = "none";
			/^-conc=(.+)$/ and $args{"concantenation"} = $1;
			/^-Gblocks=(.+)$/ and $args{"Gblocks"} = $1;
			/^-output=(.+)$/ and $args{"format"} = $1;
			/^-overwrite/ and $args{"overwrite"} = 1;
			/^-sort/ and $args{"sort"} = 1;
			/^-N/ and $args{"enize"} = 1;
			(/^-fillends/ or /^-fill/) and $args{"fillends"} = 1;
			(/^-fillgaps/ or /^-fill/) and $args{"fillgaps"} = 10;
			/^-fillgaps=(\d+)$/ and $args{"fillgaps"} = $1;
			/^-trim/ and $args{"deEmpty"} = 1;
			/^-unalign/ and $args{"unalign"} = 1;
			/^-codon/ and $args{"codon"} = 1;
			/^-numerize/ and $args{"numerize"} = 1;
			/^-translate=(.+)$/ and $args{"translate"} = $1;
			/^-partition=(.+)$/ and $args{"partition"} = $1;
			# more args to be added
		}else{
			push @files, $_;
		}
	}
}


## Verify arguments ##

die "Error: no input file specified.\n" unless @files;
foreach (@files){
	die "Error: file $_ does not exist.\n" unless -e($_);
}
die "Error: filter list ".$args{'filter'}." does not exist.\n" if ($args{"filter"} and not (-e $args{"filter"}));
die "Error: Gblocks program ".$args{'Gblocks'}." does not exist.\n" if ($args{"Gblocks"} and not (-e $args{"Gblocks"}));


#####################
#### Filter taxa ####
#####################

## Read filter list ##

my %filters;
if ($args{"filter"}){
	my $isFasta = 0;
	open FIN, "<".$args{"filter"};
	while (<FIN>){
		s/\s+$//g;
		next unless $_;
		next if /^#/;
		if (/^>/){
			$isFasta = 1;
			$_ =~ s/^>(\S+).*$/$1/;
			$filters{$_} = $_;
		}else{
			next if $isFasta;
			@a = split (/\t/);
			$filters{$a[0]} = $a[0];
		}
	}
	close FIN;
}

#############################################
#### Modification of sequence alignments ####
#############################################

foreach my $file (@files){

	my $stem = $file;
	$file =~ /(.+)\.[^.]+$/ and $stem = $1;

	## Read sequences ##
	
	my @taxa = ();
	my @seqs = ();
	my $isReading = 0;											# If this taxon should be read
	
	open FIN, "<$file";
	while (<FIN>){
		s/\s+$//g;												# chomp
		next unless $_;										# skip empty line
		next if /^#/;											# skip annotation
		if (/^>/){												# taxon name
			$_ =~ /^>(\S+).*$/;									# read until first white space
			if ($args{"filter"}){								# if the sequences should be filtered
				if (exists $filters{$1}){
					$isReading = 1;
					push (@taxa, $1);							# create new taxon in the memory
					push (@seqs, "");
				}else{
					$isReading = 0;
				}
			}else{
				$isReading = 1;
				push (@taxa, $1);
				push (@seqs, "");
			}
		}else{
			next unless $isReading;
			$_ =~ s/\s//g;										# remove white spaces
			$seqs[$#seqs] .= uc($_);
		}
	}
	close FIN;


	## Modify sequences ##

	if ($args{"deEmpty"}){										# remove empty sites (columns)
		for ($i = 0; $i <= length($seqs[0]); $i++){
			$s = 0;
			for ($j = 0; $j <= $#taxa; $j++){
				if ((substr($seqs[$j],$i,1) ne "-") && (substr($seqs[$j],$i,1) ne "N") && (substr($seqs[$j],$i,1) ne "?")){
					$s = 1;
					next;
				}
			}
			if ($s == 0){
				for ($j = 0; $j <= $#taxa; $j++){
					$seqs[$j] = substr($seqs[$j],0,$i).substr($seqs[$j],$i+1);
				}
				$i--;
			}
		}
	}

	for ($i = 0; $i <= $#taxa; $i++){
		if ($args{"unalign"}){									# remove all gaps (unalign)
			$seqs[$i] =~ s/-//g;
		}else{
			if ($args{"enize"}){								# replace ambiguous codes with 'N's
				$seqs[$i] =~ s/[RYMKWSBDHV]/N/g;
			}
			if ($args{"fillends"}){								# replace initial and final gaps with 'N's
				$seqs[$i] =~ /(^-*)([^-]+.*[^-]+)(-*$)/;
				$seqs[$i] = "N"x length($1) .$2. "N"x length($3);
			}
			if ($args{"fillgaps"}){								# replace in-sequence gaps with 'N's
				$j = $args{"fillgaps"};
				while ($seqs[$i] =~ /(^.*[^-])(-{$j,})([^-].*$)/){
					$seqs[$i] = $1."N"x length($2).$3;
				}
			}
		}
	}

	if ($args{"numerize"}){										# numerize taxon names
		open FOUT, ">$stem.translate.txt";
		for ($i = 0; $i <= $#taxa; $i++){
			$j = $i+1;
			print FOUT "taxon$j\t$taxa[$i]\n";
			$taxa[$i] = "taxon$j";
		}
		close FOUT;
	}

	if ($args{"translate"} and -e $args{"translate"}){			# translate taxon names
		my %tbTrans = ();
		open IN, "<".$args{"translate"};
		while (<IN>){
			s/\s+$//; next unless $_; next if /^#/;
			@a = split (/\t/);
			$tbTrans{$a[0]} = $a[1] if $#a;
		}
		close IN;
		for ($i = 0; $i <= $#taxa; $i++){
			if (exists $tbTrans{$taxa[$i]}){
				$taxa[$i] = $tbTrans{$taxa[$i]}
			}
		}
	}

	for ($i = 0; $i <= $#taxa; $i++){							# merge taxon names and sequences
		$taxa[$i] .= "\n".$seqs[$i]."\n";
	}
	
	if ($args{"sort"} and not $args{"numerize"}){				# sort taxon names
		@taxa = sort(@taxa);
	}

	open FOUT, ">$file.1";										# write output
	for ($i = 0; $i <= $#taxa; $i++){
		print FOUT ">$taxa[$i]";
	}
	close FOUT;

	if ($args{"Gblocks"}){										# perform Gblocks
		$s = $args{"Gblocks"}." $file.1 -t=";
		if ($args{"type"} eq "codon"){
			$s.= "c -b2=";
		}else{
			$s.= "d -b2=";
		}
		$s .= (int(($#taxa+1)/2)+1);
		$s .= " -b3=3 -b4=6 -b5=a -e=.2";
		system $s;
		open FIN, "<$file.1.2";
		open FOUT, ">$file.1";
		while (<FIN>){
			$_ =~ s/\s+$//g;
			$_ =~ s/\s//g unless (/^>/);
			print FOUT "$_\n";
		}
		close FIN;
		close FOUT;
		unlink "$file.1.2";
	}


	## Convert format ##

	if ($args{"format"} ne "fasta" and not $args{"concantenation"}){

		@taxa = (); @seqs = ();
		open FIN, "<$file.1";
		while (<FIN>){
			s/\s+$//g; next unless $_; next if /^#/;
			if (/^>(.+)$/){
				push (@taxa, $1);
				push (@seqs, "");
			}else{
				$seqs[$#seqs] .= $_;
			}
		}
		close FIN;

		if ($args{"format"} eq "nexus"){
			open FOUT, ">$stem.nex";
			print FOUT "#NEXUS\nbegin data;\n";
			print FOUT "\tdimensions ntax=". ($#taxa+1) ." nchar=". length($seqs[0]). ";\n";
			print FOUT "\tformat datatype=dna missing=? gap=-;\nmatrix\n";
			for ($i = 0; $i <= $#taxa; $i++){
				$seqs[$i] =~ s/N/\?/g;
				print FOUT "[".($i+1)."] ".$taxa[$i]."\t".$seqs[$i]."\n";
			}
			print FOUT ";\nend;\n";
			close FOUT;
			print "$file saved as $stem.nex.\n";
		}elsif ($args{"format"} eq "phylip"){
			open FOUT, ">$stem.phy";
			print FOUT " ". ($#taxa+1) ." ". length($seqs[0]) ."\n";
			for ($i = 0; $i <= $#taxa; $i++){
				if (length($taxa[$i]) >= 10){
					$s = substr($taxa[$i], 0, 10);
				}else{
					$s = $taxa[$i] . " " x (10-length($taxa[$i]));
				}
				print FOUT $s ."   ". $seqs[$i]."\n";
			}
			close FOUT;
			print "$file saved as $stem.phy.\n";
		}elsif ($args{"format"} eq "pir"){
			open FOUT, ">$stem.pir";
			for ($i = 0; $i <= $#taxa; $i++){
				print FOUT ">DL; ".$taxa[$i]."\n".$taxa[$i].".\n";
				print FOUT $seqs[$i]."*\n";
			}
			close FOUT;
			print "$file saved as $stem.pir.\n";
		}
	}


	## Write codons ##

	if ($args{"codon"}){
		@taxa = (); @seqs = ();
		open FIN, "<$file.1";
		while (<FIN>){
			s/\s+$//g; next unless $_; next if /^#/;
			if (/^>(.+)$/){
				push (@taxa, $1);
				push (@seqs, "");
			}else{
				$seqs[$#seqs] .= $_;
			}
		}
		close FIN;
		open FOUT1, ">$stem.codon1.fasta";
		open FOUT2, ">$stem.codon2.fasta";
		open FOUT3, ">$stem.codon3.fasta";
		$s = length($seqs[0]);
		$t = int($s/3);
		for ($i = 0; $i <= $#taxa; $i++){
			print FOUT1 ">$taxa[$i]\n";
			print FOUT2 ">$taxa[$i]\n";
			print FOUT3 ">$taxa[$i]\n";
			for ($j = 0; $j <= $t; $j++){
				print FOUT1 substr($seqs[$i],$j*3,1) if ($j*3<$s);
				print FOUT2 substr($seqs[$i],$j*3+1,1) if ($j*3+1<$s);
				print FOUT3 substr($seqs[$i],$j*3+2,1) if ($j*3+2<$s);
			}
			print FOUT1 "\n";
			print FOUT2 "\n";
			print FOUT3 "\n";
		}
		close FOUT1;
		close FOUT2;
		close FOUT3;
	}


	## Divide master alignment by partition ##

	if ($args{"partition"} and -e $args{"partition"}){			# translate taxon names
		@taxa = (); @seqs = ();
		open FIN, "<$file.1";
		while (<FIN>){
			s/\s+$//g; next unless $_; next if /^#/;
			if (/^>(.+)$/){
				push (@taxa, $1);
				push (@seqs, "");
			}else{
				$seqs[$#seqs] .= $_;
			}
		}
		close FIN;
		my %parts = ();
		open IN, "<".$args{"partition"};
		while (<IN>){
			s/\s+$//; next unless $_; next if /^#/;
			s/^\S+,\s//;
			@a = split (/ = /);
			$parts{$a[0]} = $a[1];
		}
		close IN;
		foreach my $part (keys %parts){
			@a = split (/, /, $parts{$part});
			my @b = ();
			foreach (@a){
				$s = 1;
				$s = $1 if s/\\(\d+)$//;
				/^(\d+)-(\d+)$/;
				for ($t=$1; $t<=$2; $t+=$s){
					push (@b, $t);
				}
			}
			open FOUT, ">$part.fasta";
			for ($i = 0; $i <= $#taxa; $i++){
				$s = "";
				foreach (@b){
					$s .= substr($seqs[$i], $_-1, 1);
				}
				print FOUT ">$taxa[$i]\n$s\n";
			}
			close FOUT;
		}
	}
}

for ($i = 0; $i <= $#files; $i++){								# overwrite original files or redefine output files
	print "$files[$i] processed";
	if ($args{"overwrite"}){
		unlink $files[$i];
		rename "$files[$i].1", $files[$i];
	}elsif ($args{"format"} ne "fasta" and not $args{"concantenation"}){
		unlink "$files[$i].1";
	}else{
		$files[$i] =~ /(.+)\.([^.]+)$/;
		$s = "$1.out.$2";
		rename "$files[$i].1", $s;
		$files[$i] = $s;
		print " and saved as $files[$i]";
	}
	print ".\n";
}


exit 0 unless $args{"concantenation"};

###########################################
## Concantenation of sequence alignments ##
###########################################

# read partition sizes

foreach my $file (@files){
	open FIN, "<$file";
	$s = $file;
	$s=~ s/\..*//; # remove extension
	$i = 0;
	while (<FIN>){
		s/\s+$//g;
		next unless $_;
		last if exists $parts{$s};
		if (/^>/){
			if ($i){
				push @part_order, $s;
				$parts{$s} = $i;
			}
		}else{
			s/\s//g;
			$i += length($_);
		}
	}
	close FIN;
}
print (($#files+1)." partitions read.\n");


## Read taxon names ##

foreach my $file (@files){
	open FIN, "<$file";
	while (<FIN>){
		s/\s+$//g;
		/^>/ or next;
		s/^>(\S+).*$/$1/;
		$seqs{$_} = "";
	}
	close FIN;
}
print keys(%seqs) ." taxa read.\n";


## Read sequences ##

$i = 0; # accumulated sequence length
foreach my $file(@files){
	open FIN, "<$file";
	while (<FIN>){
		s/\s+$//g;
		next unless $_;
		if (/^>/) {
			# $_ =~ s/^>(\S+).*$/$1/; # get the sequence name #### THIS IS MY TOP SECRET.
			/^>(\S+)/;
			$s = $1;
		}else{
			if (exists $seqs{$s}){
				s/\s//g; # remove spaces;
				$seqs{$s} .= $_;
			}
		}
	}
	close FIN;
	$s = $file;
	$s=~ s/\..*//;
	$i += $parts{$s};
	if ($args{"fillends"}){
		$t = "N"x$parts{$s};
	}else{
		$t = "-"x$parts{$s};
	}
	foreach $key (keys %seqs){
		if (length($seqs{$key}) < $i){
			$seqs{$key} .= $t;
		}
	}
	print "Partition $s processed.\n";
	unlink $file;
}

## Output in different formats ##

if ($args{"concantenation"} eq "none"){										# default (fasta)
	open FOUT, ">output.fasta";
	foreach $key (sort (keys %seqs)){
		print FOUT ">$key\n$seqs{$key}\n";
	}
	close FOUT;
	print "The concantenated alignment has been saved as output.fasta.\n";

}elsif ($args{"concantenation"} eq "raxml"){									# RAxML
	open FOUT, ">output.phy";
	print FOUT " ". keys (%seqs) ." ";
	foreach $key (keys %seqs){
		print FOUT length($seqs{$key})."\n";
		last;
	}
	foreach $key (sort (keys %seqs)){
		if (length($key) >= 10){
			$s = substr($key, 0, 10);
		}else{
			$s = $key . " " x (10-length($key));
		}
		print FOUT $s ."   ". $seqs{$key}."\n";
	}
	close FOUT;
	print "The concantenated alignment has been saved as output.phy.\n";
	if ($#files){ # write partitioning information
		open FOUT, ">output_partitions.txt";
		$i = 0;
		foreach $key (@part_order){
			$j = $i + $parts{$key};
			if ($key =~ /^\d\dS$/){ # ribosomal RNA genes
				print FOUT "DNA, $key = ". ($i+1) ."-$j\n";
			}else{ # coding genes
				print FOUT "DNA, $key-1 = ". ($i+1) ."-$j\\3\n";
				print FOUT "DNA, $key-2 = ". ($i+2) ."-$j\\3\n";
				print FOUT "DNA, $key-3 = ". ($i+3) ."-$j\\3\n";
			}
			$i = $j;
		}
		close FOUT;
		print "The partitioning information has been saved as output_partitions.txt.\n";
	}
	
}elsif ($args{"concantenation"} eq "beast"){									# BEAST
	open FOUT, ">output.nex";
	print FOUT "#NEXUS\nbegin taxa;\n";
	print FOUT "\tdimensions ntax=". keys(%seqs) .";\n\ttaxlabels\n";
	foreach $key (sort (keys %seqs)){
		print FOUT "\t$key\n";
	}
	print FOUT ";\nend;\n\n";
	print FOUT "begin characters;\n";
	print FOUT "\tdimensions nchar=". $i .";\n";
	print FOUT "\tformat datatype=dna missing=? gap=-;\n\tmatrix\n";
	foreach $key (sort (keys %seqs)){
		print FOUT "\t$key\t";
		$seqs{$key} =~ s/N/\?/g;
		print FOUT $seqs{$key} . "\n";
	}
	print FOUT ";\nend;\n\n";
	# partitioning information
	$i = 0;
	print FOUT "begin assumptions;\n";
	foreach $key (@part_order){
		print FOUT "\tcharset $key = ". ($i+1) ."-". ($i+=$parts{$key}) .";\n";
	}
	print FOUT "end;\n";
	close FOUT;
	print "The concantenated alignment has been saved as output.nex.\n";

}elsif ($args{"concantenation"} eq "mrbayes"){									# MrBayes
	open FOUT, ">output.nex";
	print FOUT "#NEXUS\nbegin data;\n";
	print FOUT "\tdimensions ntax=". keys(%seqs) ." nchar=". $i .";\n";
	print FOUT "\tformat datatype=dna missing=? gap=-;\n\tmatrix\n";
	foreach $key (sort (keys %seqs)){
		print FOUT "\t$key\t";
		$seqs{$key} =~ s/N/\?/g;
		print FOUT $seqs{$key} . "\n";
	}
	print FOUT ";\nend;\n\n";
	$i = 0; $s = "";
	print FOUT "begin mrbayes;\n";
	foreach $key (@part_order){
		$s .= ", " if $s;
		$s .= $key;
		print FOUT "\tcharset $key = ". ($i+1) ."-". ($i+=$parts{$key}) .";\n";
	}
	print FOUT "\tpartition scheme1=".@part_order.": $s;\n";
	print FOUT "\tset partition=scheme1;\n";
	print FOUT "\tlset app=(all) nst=mixed rates=invgamma;\n";
	print FOUT "\tunlink revmat=(all) pinvar=(all) shape=(all) statefreq=(all);\n";
	print FOUT "\tprset ratepr=variable;\n";
	print FOUT "end;\n";
	close FOUT;
	print "The concantenated alignment has been saved as output.nex.\n";

}

print "Task completed.\n";
exit 0;


