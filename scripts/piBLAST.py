#!/usr/bin/env python
# this script checks if any piRNA sequences BLAST to the TE sequences with varying percent identities
# first 8 bases of piRNA must match
# USE: piBLAST.py

import re
import sys
import os
from subprocess import Popen, PIPE 
from collections import defaultdict
from collections import Counter
import pickle
import itertools

pi_IN="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/WB_piRNA_positions.gff"
reference="/lscr2/andersenlab/kml436/sv_sim2/c_elegans.PRJNA13758.WS245.genomic.fa"
pi_fasta="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/piRNAs.fasta"
TE_consensus="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/SET2/round2_consensus_set2.fasta"
family_renames="/lscr2/andersenlab/kml436/git_repos2/Transposons2/files/round2_WB_familes_set2.txt"


# put shortened WB family names into a dictionary
renames={}
with open(family_renames, 'r') as IN:
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('\t',line)
		element,family=items[0:2]
		renames[element]=family

# make blast database of TE sequences if it doesn't already exist
if not os.path.isfile("TE_database.nsq"):
	cmd="/lscr2/andersenlab/kml436/ncbi-blast-2.2.30+/bin/makeblastdb -in  {TE_consensus} -dbtype nucl -out TE_database".format(**locals())
	result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()
else:
	print "BLAST database already exists, continuing..."

# blast piRNA seqeunces to TE sequences
cmd="/lscr2/andersenlab/kml436/ncbi-blast-2.2.30+/bin/blastn -db TE_database -query {pi_fasta} -evalue .1 -word_size 5 -outfmt '6 qseqid sseqid pident qlen length mismatch gapopen evalue bitscore qstart qend btop' -max_target_seqs 100 -out piRNA_blast.txt -num_threads 10".format(**locals())
result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()



OUT=open("piRNA_blast_strict_redundant.txt", 'w')
with open("piRNA_blast.txt".format(**locals()) ,'r') as IN:
	for line in IN:
		line=line.rstrip()
		items=re.split('\t',line)
		TE=items[1]
		query_start=int(items[9])
		btop=items[11]
		btop_nums=re.findall('\d+', btop)
		first_digit=int(btop_nums[0])
		if TE in renames.keys():
			TE=renames[TE]
			items[1]=TE
		new_line='\t'.join(items[0:])
		#if query_start==1 : #and first_digit>=8
		OUT.write(new_line + '\n')
OUT.close()

cmd="cat piRNA_blast_strict_redundant.txt |sort -k1,1 -k2,2 -k10,10 -k11,11r -k8,8 -k9,9 -k3,3r |uniq  > piRNA_blast_strict_21.txt" #| awk '$11>20 {print $0}'
result, err = Popen([cmd],stdout=PIPE, stderr=PIPE, shell=True).communicate()


seen={}

OUT=open("summary_mismatches_BLAST_strict.txt", 'w')
OUT.write("Number of Mismatches\tNumber Unique piRNAs Aligned to One TE\tNumber Unique piRNAs Aligned to Multiple TEs\n")

BLAST_PAIRS=open("blast_pairs.txt", 'w')

mis_per= {'zero': 100, 'one': 95.23,'two': 90.48, 'three': 85.71,'four':80.95,'five':76.19}
num_ver= {'zero': 0, 'one': 1,'two': 2, 'three': 3, 'four':4, 'five':5}
def piblast(mismatch):
	blasts={}
	blasts=defaultdict(list)

	pi_one=0
	pi_multiple=0
	with open("piRNA_blast_strict_21.txt", 'r') as IN:
		for line in IN:
			print line
			line=line.rstrip('\n')
			items=re.split('\t',line)
			query,TE,perID=items[0:3]
			match = re.search("(?:Pseudogene|Transcript|sequence_name|^Name)(?:=|:)([\w|\d]+.\d+)", query) #just pull gene name, remove splice info
			pi_transcript =match.group(1)
			perID=items[2]
			beat_per=mis_per[mismatch]


			length_align=items[10]
			info_align=items[11]
			matches=re.split("\D+",info_align)
			matches=[int(i) for i in matches]
			matched_bases=sum(matches)
			print matched_bases



			num=num_ver[mismatch]
			print num
			family_short=re.sub("_CE$","",TE)
			family_short=re.sub("WBTransposon","WBT",family_short)
			pair=family_short + "_" + pi_transcript
			actual_mismatch=21-int(matched_bases)
			print actual_mismatch
			if actual_mismatch<=num:
			#if float(perID)>=beat_per: 
				if pair not in seen.keys():
					blasts[family_short].append(pi_transcript)
					BLAST_PAIRS.write("{pi_transcript}\t{family_short}\t{num}\n".format(**locals()))
					seen[pair]=mismatch


	blasts_TEs_strict = len(blasts.keys())
	vals=list(itertools.chain(blasts.values()))
	blasts_pis_strict =len(set(list(itertools.chain.from_iterable(vals))))

	#ditct[TE].pi,pi,pi
	blasts_pis_new=list(itertools.chain.from_iterable(vals))
	pi_counts=Counter(blasts_pis_new)
	print pi_counts


	for k,v in pi_counts.items():
		if v==1:
			pi_one+=1
		else:
			pi_multiple+=1

	#OUT.write("{mismatch}\t{blasts_pis_strict}\t{blasts_TEs_strict}\n".format(**locals()))
	OUT.write("{mismatch}\t{pi_one}\t{pi_multiple}\n".format(**locals()))


	with open("strict_blasts_{mismatch}.txt".format(**locals()), "wb") as fp: # Pickle
		pickle.dump(blasts, fp) 
	


piblast('zero')
piblast('one')
piblast('two')
piblast('three')
piblast('four')
piblast('five')
#for k,v in seen.items():
#	print k
#	print v
OUT.close()
BLAST_PAIRS.close()

