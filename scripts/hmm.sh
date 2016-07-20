#!/bin/bash
files_dir=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files
#CDS
cds_arms=`bedtools intersect -wao -a ${files_dir}/arms.gff -b ${files_dir}/WB_CDS_positions.gff |datamash sum 19`
cds_centers=`bedtools intersect -wao -a ${files_dir}/centers.gff -b ${files_dir}/WB_CDS_positions.gff |datamash sum 19`
#Promoter
promoter_arms=`bedtools intersect -wao -a ${files_dir}/arms.gff -b ${files_dir}/WB_promoter_positions.gff |datamash sum 19`
promoter_centers=`bedtools intersect -wao -a ${files_dir}/centers.gff -b ${files_dir}/WB_promoter_positions.gff |datamash sum 19`

intron_arms=`bedtools intersect -wao -a ${files_dir}/arms.gff -b ${files_dir}/WB_intron_positions.gff |datamash sum 19`
intron_centers=`bedtools intersect -wao -a ${files_dir}/centers.gff -b ${files_dir}/WB_intron_positions.gff |datamash sum 19`

gene_arms=`bedtools intersect -wao -a ${files_dir}/arms.gff -b ${files_dir}/WB_gene_positions.gff |datamash sum 19`
gene_centers=`bedtools intersect -wao -a ${files_dir}/centers.gff -b ${files_dir}/WB_gene_positions.gff |datamash sum 19`

echo cds_arms "<-" $cds_arms
echo cds_centers "<-" $cds_centers
echo cds_total "<-" `expr $cds_arms + $cds_centers`
echo promoter_arms "<-" $promoter_arms
echo promoter_centers "<-" $promoter_centers
echo promoter_total "<-" `expr $promoter_arms + $promoter_centers`
echo intron_arms "<-" $intron_arms
echo intron_centers "<-" $intron_centers
echo intron_total "<-" `expr $intron_arms + $intron_centers`
#echo gene arms: $gene_arms
#echo gene centers: $gene_centers

intergenic_arms=`expr $gene_arms - $promoter_arms`
intergenic_centers=`expr $gene_centers - $promoter_centers`

echo intergenic_arms "<-" $intergenic_arms
echo intergenic_centers "<-" $intergenic_centers
echo intergenic_total "<-" `expr $intergenic_arms + $intergenic_centers`





#all[7,]<-c("Intergenic",ctest$statistic,ctest$p.value)


