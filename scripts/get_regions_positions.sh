#!/bin/bash
files_dir=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/

cds_arms=`bedtools intersect -wao -a ${files_dir}/arms.gff -b ${files_dir}/WB_CDS_positions.gff |datamash sum 19`

echo $cds_arms