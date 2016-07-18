#!/usr/bin/env python
# this script calculates cumulative totals of insertion and reference calls for all traits
# USE: overall_te.py <T_kin_C_matrix_full>

import sys
import re
from collections import defaultdict

insertion_traits=defaultdict(list)
reference_traits=defaultdict(list)
cumulative=defaultdict(list)

OUT=open("per_strain_NAs_removed.txt",'w')

infile=sys.argv[1]

with open(infile, 'r') as IN:
	header=next(IN)
	for line in IN:
		line=line.rstrip('\n')
		items=re.split('\t',line)
		trait= items[0]
		if trait !="coverage":
			id_info=re.split("_TRANS_",trait)
			
			method=id_info[0]
			TE=id_info[1]


			if method=="new" or method== "ONE_new":
				strain_info= items[1:]
				insertion_traits[TE]=strain_info
			elif method=="reference":
				strain_info= items[1:]
				reference_traits[TE]=strain_info
				#print id_info
				#print method
				#print TE




common_TEs = set(insertion_traits.keys()) & set(reference_traits.keys())


for i in common_TEs:
	for strain in range(0,len(insertion_traits[i])):
		ins=(insertion_traits[i])[strain]
		ref=(reference_traits[i])[strain]
		if ins !="NA" and ref!="NA":

			total=float(ins)+float(ref)

			cumulative[i].append(total)
		else:
			cumulative[i].append("NA")


#cumulative[key] = map("cumulative_TRANS"+cumulative[key], cumulative[key])

cumulative={"cumulative_TRANS_" + key:value for key,value in cumulative.items()}

print cumulative.keys()

OUT.write("trait")
for i in (re.split('\t',header))[1:]:
	OUT.write('\t' + i)
OUT.write('\n')

for i,v in cumulative.items():
	OUT.write(i)
	for strain in v:
		OUT.write('\t' +str(strain))
	OUT.write('\n')

#result=[ sum(x) for x in zip(insertion_traits[i], reference_traits[i])] 

OUT.close()

sys.exit()

uniq_to_ins = set(insertion_traits.keys()) - set(reference_traits.keys())
uniq_to_ref = set(reference_traits.keys()) - set(insertion_traits.keys())

print len(insertion_traits)
print len(reference_traits)

print len(common_TEs)
print common_TEs
print '\n'

print len(uniq_to_ins)
print uniq_to_ins
print '\n'

print len(uniq_to_ref)
print uniq_to_ref
print '\n'
