#!/usr/bin/env python

from subprocess import Popen, PIPE
import re
from collections import defaultdict
gene=defaultdict(list)

gene_gff="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/WB_gene_positions.gff"
taj="/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/TajimaD_interest.bed"


result, err = Popen(["""bedtools intersect -wao -a {taj} -b {gene_gff} |cut -f1,2,10,11,15> /lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/pre_taj_genes.txt""".format(**locals())], stdout=PIPE, stderr=PIPE, shell=True).communicate()


OUT=open("/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/taj_genes.txt",'w')
#OUT.write("Chromosome\tBin Start\tGene Start\tGene End\tGene\n")
with open("/lscr2/andersenlab/kml436/git_repos2/Transposons2/results/final_results/pre_taj_genes.txt", 'r') as IN:
	for line in IN:
		line=line.rstrip()
		items=re.split('\t',line)
		chrom,bin_s,gene_s,gene_e,gene_info=items[0:5]
		gene_info=re.split(';',gene_info)
		ID=chrom+bin_s
		seq_name=re.sub("sequence_name=","",gene_info[3])

		gene[ID].append(seq_name)
		#OUT.write("{chrom}\t{bin_s}\t{gene_s}\t{gene_e}\t{seq_name}\n".format(**locals()))
for k,v in gene.items():
	genes=','.join(v[0:len(v)])
	OUT.write("{k}\t{genes}\n".format(**locals()))



		#OUT.write("{chrom}\t{bin_s}\t{gene_s}\t{gene_e}\t{seq_name}\n".format(**locals()))

OUT.close()

