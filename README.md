BeforePhylo
==========

a Perl script for manipulating multiple sequence alignments for phylogenetic reconstruction

Version 0.9.0

by Qiyun Zhu (<qiyunzhu@gmail.com>)

## Features

* Concatenate multiple MSAs and generate partition table for RAxML, MrBayes and BEAST
* Remove empty sites, fill gaps and/or ends, replace ambiguous codes with 'N's
* Translate taxon names into numbers or based on a dictionary
* Extract subset of sequences based on a list
* Divide MSA by codon position or by a partition table
* Convert formats (FASTA, NEXUS and Phylip)

## Usage

    perl BeforePhylo.pl [option(s)] [MSA(s)]

### Input:

One or more MSA files in FASTA format.

### Options:

- `-type=<dna|protein|condon>`

    Type of data (default=dna).

- `filter=<taxa list file>`

    Only retain sequences defined in a taxa list file.

- `conc=<none|raxml|beast|mrbayes>`

    Concantenate multiple MSAs and generate partition table.

- `sort`

    Sort sequence names in alphabetical order.

- `Gblocks=<Gblocks program>`

    Perform Gblocks on each MSA.

- `trim`

    Remove empty sites.

- `fillends`

    Replace initial and final gaps with 'N's.

- `fillgaps=<cutoff, default=10>`

    Replace in-sequence gaps no shorter than the cutoff with 'N's.

- `fill`

    Fill both ends and gaps.

- `unalign`

    Remove all gaps (make the sequences unaligned).

- `N`

    Replace ambiguous codes with 'N's.

- `10`

    Keep the first 10 characters of sequence names (for old-style Phylip format).

- `codon`

    Divide dataset into three codon positions.

- `numerize`

    Translate taxon names into numbers, to avoid too complicated / duplicated names.

- `translate=<dictionary file>`

    Translate taxon names according to a dictionary.

- `partition=<RAxML-style partition table>`

    Split a master alignment by partition.

- `output=<fasta|nexus|phylip|pir>

    Output file format (default=fasta).

- `overwrite`

    Overwrite original files (default is to create new files).


## Examples

  Concatenate three gene markers into one master alignment, and create a RAxML-style partition table:
  
    perl BeforePhylo.pl -output=phylip -conc=raxml coI.fas coII.fas cytB.fas

  Clean up all alignments in the current directory, removing empty sites, filling terminal spaces, and sort sequences by name:
  
    perl BeforePhylo.pl -N -trim -fillends -sort -overwrite *.fas

  For numerous, long, duplicated and space-containing sequences names (depreciated by many programs), convert them into numbers. After tree-building, you can use AfterPhylo.pl to convert them back:
  
    perl BeforePhylo.pl -numerize 16S.fas

  Create three alignments, eaching containing one codon of the original alignment. They are subject to model test and/or phylogenetic reconstruction.

    perl BeforePhylo.pl -codon groEL.fas

  Extract a subset of taxa based on a user-defined list:

    perl BeforePhylo.pl -filter=taxa_list.txt mtDNA.fas

  The file `taxa_list.txt` looks like:
  
    Dmel
    Dsim
    Dpse
    Dvir
    ...


## Dependencies

This program does not depend on any third-party modules.


## License

Licensed under the [BSD 2-clause license](http://opensource.org/licenses/BSD-2-Clause).
