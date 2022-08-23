# yexinhai, yexinhai@zju.edu.cn
#Run orthofinder for orthology inference.
orthofinder -f ./01_data -t 112 -a 10 -M msa -T fasttree
# Generate a supergene sequence from 1,792 single-copy orthologues.
python concat.py
#Run IQ-TREE v2.0 for maximum likelihood phylogenetic tree construction.
iqtree2 -s ./all_concated.fasta -m MFP -B 1000
#Run r8s
r8s -b -f r8s_input.txt >r8s.out
#Run CAFE v5.0
cafe5 -i gene_families.txt -t r8s.tre



#Run EDTA pipeline for TEs identification.
perl EDTA.pl -genome $genome -step all -sensitive 1 -anno 1 -t 112
python time.py # Insertion times estimate for full-length LRT-RTs

#run insectOR pipeline for identication of OR genes
python split.py
python run_exonerate.py
python merge.py
perl scoreGenesOnScaffold.pl -i exonerate.txt -s $genome -q query.fasta -tmh1 -tmh2 -tmh3 -p -m &
mafft --maxiterate 1000 --localpair --thread 112 --quiet all.ComOR.pep.fa >all.ComOR.pep.mafft
trimal -in all.ComOR.pep.mafft -out all.ComOR.pep.trimal -automated1
iqtree2 -s all.ComOR.pep.trimal -m MFP -B 1000 -T AUTO
