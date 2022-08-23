#!/usr/bin/python
# yexinhai, yexinhai@zju.edu.cn

import os
import pandas as pd
import re
from Bio import SeqIO

filepath = '/data/yangyi/project/16_Anastatus/03_genomic_feature/01_EDTA'
species = ['Acep','Amel','Aros','Bter','Btre','Cchi','Cflo','Csol','Dcol','Gfla','Nvit','Obir']

def split2file(LTRs) :
	adict = {}
	for line in open(LTRs,'r') :
		line = line.strip()
		if line.startswith('>') :
			seqid = line.split(' ')[1]
			adict[seqid] = ''
		else :
			adict[seqid] += line
	Ls = list(set([i[1:] for i in list(adict.keys())]))
	for i in Ls:
		w = open('./tmp/' + i + '.fa','a')
		w.write('>%s\n%s\n' %('l'+ i, adict['l' + i]))
		w.write('>%s\n%s\n' %('r' + i , adict['r' + i]))
		w.close()
		os.system('mafft --quiet --thread 112 --maxiterate 1000 --localpair ./tmp/%s.fa >./tmp/%s.mafft' %(i,i))
		os.system('distmat -nucmethod 2 ./tmp/%s.mafft -outfile ./tmp/%s.distmat' %(i,i))

for s in species :
	os.system('mkdir %s' %s)
	os.chdir('./%s' %s)
	os.system('ln -s %s/%s.genome.fa.mod %s.genome.fa' %(filepath,s,s))
	os.system('ln -s %s/%s.genome.fa.mod.EDTA.intact.gff3 ./' %(filepath,s))
	os.system('gfftobed -f long_terminal_repeat -a ID %s.genome.fa.mod.EDTA.intact.gff3 >LTRs.bed' %s)
	os.system('sed -i "s/?/+/" LTRs.bed')
	os.system('seqkit subseq --bed  LTRs.bed -o LTRs.fa --id-ncbi  %s.genome.fa' %s)
	os.system('mkdir tmp')
	split2file('LTRs.fa')

	outlist = []
	for file in os.listdir('./tmp') :
		if file.endswith('.distmat') :
			with open('./tmp/' + file ,'r') as f :
				lines = f.readlines()
			outlist.append([file.split('.')[0],float(lines[8].split('\t')[2].strip())])
	
	df = pd.DataFrame(outlist,columns=['id','d'])
	
	def get_te_name(d,gff) :
		adict = {}
		for line in open(gff,'r') :
			line = line.strip()
			a = line.split('\t')
			if a[2] == 'long_terminal_repeat' :
				i = re.search('ID=(\S+?);',line).group(1)[1:]
				n = re.search('Name=(\S+?);',line).group(1)
				c = re.search('Classification=(\S+?);',line).group(1)
				adict[i] = [n,c]
		return(adict[d[0]])
	df[['Name','Class']] = df.apply(lambda x : get_te_name(x,'./%s.genome.fa.mod.EDTA.intact.gff3' %s),axis=1,result_type='expand')
	df = df[['Name','id','Class','d']]
	df.to_excel('%s.intactLTR.xlsx' %s,index=False)
	os.chdir('../')
