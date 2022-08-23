#!/usr/bin/python
# yexinhai, yexinhai@zju.edu.cn

import pandas as pd
from Bio import SeqIO
import os
from collections import defaultdict


seqindex = SeqIO.index('./01_data/all.pep.fa','fasta')

og2genes = {}
with open('Orthogroups.txt','r') as f :
	lines = f.readlines()
for line in lines :
	line = line.strip()
	og = line.split(':')[0]
	genes = line.split(':')[1].lstrip(' ').split(' ')
	og2genes[og] = genes

data = pd.read_csv('Orthogroups.GeneCount.tsv',header=0,sep ='\t')
data.drop('Total',axis=1,inplace=True)
c = data.columns.tolist()[1:]

for s in c :
	data = data[data[s]==1]
og111 = data['Orthogroup'].tolist()

os.system('mkdir 111og')

outdict = defaultdict(list)
for og in og111 :
	genes = og2genes[og]
	w = open('./111og/' + og + '.fa','a')
	for gene in genes :
		w.write('>%s\n%s\n' %(gene.split('|')[0],str(seqindex[gene].seq)))
	w.close()
	os.system('mafft --thread 56 --anysymbol --auto ./111og/%s >./111og/%s.mafft' %(og + '.fa',og))
	os.system('trimal -in ./111og/%s -out ./111og/%s.trimal -automated1' %(og + '.mafft',og))
	for line in open('./111og/%s.trimal' %og,'r') :
		line = line.strip()
		if line.startswith('>') :
			id = line[1:]
		else :
			outdict[id].append(line.replace('\n',''))

out = open('all_concated.fasta','a+')
for key in list(outdict.keys()) :
	out.write('>' + key)
	out.write('\n')
	out.write(''.join(outdict[key]))
	out.write('\n')
out.close()
