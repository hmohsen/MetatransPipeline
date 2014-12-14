#!/bin/bash

# $1 = sample; $2 = velvetDirectory; $3 = OasesDirectory; $4 = resultsDirectory;

mkdir $4

$2./velveth $4/res_directory_25 25 -fastq -short $1_1.fastq $1_2.fastq -strand_specific
$2./velvetg $4/res_directory_25 -read_trkg yes -amos_file yes
$3./oases $4/res_directory_25

$2./velveth $4/res_directory_27 27 -fastq -short $1_1.fastq $1_2.fastq -strand_specific
$2./velvetg $4/res_directory_27 -read_trkg yes
$3./oases $4/res_directory_27

$2./velveth $4/res_directory_29 29 -fastq -short $1_1.fastq $1_2.fastq -strand_specific
$2./velvetg $4/res_directory_29 -read_trkg yes
$3./oases $4/res_directory_29

$2./velveth $4/res_mergedAssembly 27 -long $4/res_directory*/transcripts.fa
$2./velvetg $4/res_mergedAssembly -read_trkg yes -conserveLong yes
$3./oases $4/res_mergedAssembly -merge yes
