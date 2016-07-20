#!/bin/bash
# this script analyzes GWAS results to investigate invlovement of genes of interest and piRNAs
# USE: PostMappings.sh
# NOTE: must first transfer "Peak_Table.txt" and "TOI.txt" to final_results directory

scripts_dir=/lscr2/andersenlab/kml436/git_repos2/Transposons2/scripts
data_dir=/lscr2/andersenlab/kml436/git_repos2/Transposons2/data
master_list=/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/master_sample_list.txt

# check if QTL overlap regions of interest
echo "Investigating QTL overlaps..."
python ${scripts_dir}/QTL_control_overlap.py

# check if piRNA align to TE sequences
echo "Investigating piRNA alignments to TEs..."
cd /lscr2/andersenlab/kml436/git_repos2/Transposons2/piRNA
mkdir bwa
cd bwa
python ${scripts_dir}/piBWA.py

# check alignment of piRNA to respective TE families
echo "Investigating QTL specific piRNA alignments..."
python ${scripts_dir}/pi_align_check.py

# BLAST all piRNA to all TE seqs
echo "Investigating piRNA blasted to TE seqs..."
cd /lscr2/andersenlab/kml436/git_repos2/Transposons2/piRNA
mkdir blast
cd blast
python ${scripts_dir}/piBLAST.py

# check BLASTs of piRNA to respective TE families
echo "Investigating QTL specific piRNA blasts..."
cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/piRNA/bwa/TC3_TE_seqs.fasta .
python ${scripts_dir}/pi_blast_check.py

# compare BWA and BLAST results
cd /lscr2/andersenlab/kml436/git_repos2/Transposons2/piRNA
cat bwa/BWA_pairs.txt | sort | uniq > sorted_BWA_pairs.txt
cat blast/blast_pairs.txt |sort| uniq > sorted_blast_pairs.txt

comm sorted_BWA_pairs.txt  sorted_blast_pairs.txt -1 -2 > both.txt
comm sorted_BWA_pairs.txt  sorted_blast_pairs.txt -3 -2 > BWA_only.txt
comm sorted_BWA_pairs.txt  sorted_blast_pairs.txt -3 -1 > blast_only.txt

cat both.txt | awk '{print $1"\t"$2"\t"$3"\tBoth"}' > tmp && mv tmp both.txt
cat BWA_only.txt | awk '{print $1"\t"$2"\t"$3"\tBWA Only"}' > tmp && mv tmp BWA_only.txt
cat blast_only.txt | awk '{print $1"\t"$2"\t"$3"\tBLAST Only"}' > tmp && mv tmp blast_only.txt

cat both.txt BWA_only.txt blast_only.txt | sort -k3,3 -k1,1 -k2,2 -k4,4 > table_piRNAs.txt

# move tables to "tables" directory
mkdir tables/
mv table_piRNAs.txt tables/
cp bwa/summary_mismatches_BWA_strict.txt tables/
cp blast/summary_mismatches_BLAST_strict.txt tables/


# move all tables
cd /lscr2/andersenlab/kml436/git_repos2/Transposons2/
mkdir tables
cd tables

cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/piRNA/tables/summary_mismatches_BLAST_strict.txt . 
cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/piRNA/tables/table_piRNAs.txt .
cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/piRNA/tables/summary_mismatches_BWA_strict.txt .
cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/TE_seqs.txt .
cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/gene_interrupt/essentiality_nonredundant_strain_info.txt .
cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/T_kin_matrix_full.txt .
cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/T_with_monomorphic_kin_matrix_full.txt .
cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/T_kin_C_matrix_full_reduced.txt .
cp /lscr2/andersenlab/kml436/git_repos2/Transposons2/files/WB_pos_element_names.gff .


echo -e "Chromosome\tChromosome Start\tChromosome End\tTransposable Element Name\txref\tOrientation" | cat - WB_pos_element_names.gff |cut -f1-4,6 > WB_pos_table.txt
echo -e "piRNA Transcript\tTransposon\tNumber of Mismatches\tSupporting Alignment Method" | cat - table_piRNAs.txt > tmp && mv tmp table_piRNAs.txt
cat T_kin_C_matrix_full_reduced.txt |sed 's/trait/Trait/' > tmp && mv tmp T_kin_C_matrix_full_reduced.txt
cat T_with_monomorphic_kin_matrix_full.txt |sed 's/trait/Marker/' > tmp && mv tmp T_with_monomorphic_kin_matrix_full.txt
cat T_kin_matrix_full.txt |sed 's/trait/Marker/' > tmp && mv tmp T_kin_matrix_full.txt
sed -e 's/TE_start/Transposon Position/' -e 's/Transcript_Name/Transcript Name/' -e 's/Gene_Name/Gene Name/' -e 's/GO_Annotation/GO Annotation/' -e 's/No_Strains/Number of Strains/' -e 's/TE/Transposon/' essentiality_nonredundant_strain_info.txt | cut -f1-2,4-12 >tmp && mv tmp essentiality_nonredundant_strain_info.txt
cat summary_mismatches_BWA_strict.txt|sed -e 's/Number Unique piRNAs Aligned/Number Unique piRNAs Aligned with BWA/' -e 's/Number Unique Transposons/Number Unique Transposons Aligned with BWA/' >tmp1
cat summary_mismatches_BLAST_strict.txt|sed -e 's/Number Unique piRNAs BLASTED/Number Unique piRNAs Aligned with BLAST/' -e 's/Number Unique Transposons/Number Unique Transposons Aligned with BLAST/' > tmp2
paste tmp1 tmp2 |cut -f1-3,5-6 >tmp3
mv tmp3 BWA_and_BLAST_table.txt
rm tmp1
rm tmp2
rm tmp3





mkdir final_tables
cp T_kin_C_matrix_full_reduced.txt final_tables/
cp T_kin_matrix_full.txt final_tables final_tables/
cp T_with_monomorphic_kin_matrix_full.txt final_tables/
cp essentiality_nonredundant_strain_info.txt final_tables/
cp TE_seqs.txt final_tables/
cp summary_mismatches_BLAST_strict.txt final_tables/
cp summary_mismatches_BWA_strict.txt final_tables/
cp table_piRNAs.txt final_tables/
cp WB_pos_table.txt final_tables/

mkdir raw_strain_calls
cd raw_strain_calls
while read line; do 
  echo $line
  echo "copying raw output for $line"
  mkdir $line
  cp ${data_dir}/${line}/final_results/* $line
done <$master_list




