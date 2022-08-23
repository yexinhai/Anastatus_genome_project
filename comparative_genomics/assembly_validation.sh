#CCS reads
minimap2 -ax map-hifi -t 112 -o Ajap.ccs.minimap2.sam  Ajap.genome.LG.fa Ajap.ccs.fastq.gz
samtools view -F4 -bSh -@ 112 -o Ajap.ccs.minimap2.bam Ajap.ccs.minimap2.sam
samtools sort -@112 Ajap.ccs.minimap2.bam >Ajap.ccs.minimap2.sorted.bam
samtools index -@56 Ajap.ccs.minimap2.sorted.bam
samtools view -F3844 -bSh -@ 112 -o Ajap.ccs.minimap2.primaryOnly.bam Ajap.ccs.minimap2.sam
samtools sort -@112 Ajap.ccs.minimap2.primaryOnly.bam >Ajap.ccs.minimap2.primaryOnly.sorted.bam
samtools index -@56 Ajap.ccs.minimap2.primaryOnly.sorted.bam
minimap2 -ax map-hifi -t 112 -o Aful.ccs.minimap2.sam  Aful.genome.LG.fa Aful.ccs.fastq.gz
samtools view -F4 -bSh -@ 112 -o Aful.ccs.minimap2.bam Aful.ccs.minimap2.sam
samtools sort -@112 Aful.ccs.minimap2.bam >Aful.ccs.minimap2.sorted.bam
samtools index -@56 Aful.ccs.minimap2.sorted.bam
samtools view -F3844 -bSh -@ 112 -o Aful.ccs.minimap2.primaryOnly.bam Aful.ccs.minimap2.sam
samtools sort -@112 Aful.ccs.minimap2.primaryOnly.bam >Aful.ccs.minimap2.primaryOnly.sorted.bam
samtools index -@56 Aful.ccs.minimap2.primaryOnly.sorted.bam
samtools view -@ 112 Aful.ccs.minimap2.sorted.bam |cut -f 1 > Aful.ccs.minimap2.AllReadName
sort Aful.ccs.minimap2.AllReadName |uniq -iu >Aful.ccs.minimap2.UniqueReadName
grep -w -f Aful.ccs.minimap2.UniqueReadName Aful.ccs.minimap2.sam >Aful.ccs.minimap2.unique.NoHeader.sam
grep "^@" Aful.ccs.minimap2.sam > Aful.ccs.minimap2.sam.header
cat Aful.ccs.minimap2.sam.header Aful.ccs.minimap2.unique.NoHeader.sam >Aful.ccs.minimap2.unique.sam
samtools view -@ 112 Ajap.ccs.minimap2.sorted.bam |cut -f 1 > Ajap.ccs.minimap2.AllReadName
sort Ajap.ccs.minimap2.AllReadName |uniq -iu >Ajap.ccs.minimap2.UniqueReadName
grep -w -f Ajap.ccs.minimap2.UniqueReadName Ajap.ccs.minimap2.sam >Ajap.ccs.minimap2.unique.NoHeader.sam
grep "^@" Ajap.ccs.minimap2.sam > Ajap.ccs.minimap2.sam.header
cat Ajap.ccs.minimap2.sam.header Ajap.ccs.minimap2.unique.NoHeader.sam >Ajap.ccs.minimap2.unique.sam
#ONT reads
minimap2 -ax map-ont -t 112 -o Ajap.ont.minimap2.sam  Ajap.genome.LG.fa Ajap.ont.fastq.gz
samtools view -F4 -bSh -@ 112 -o Ajap.ont.minimap2.bam Ajap.ont.minimap2.sam
samtools sort -@112 Ajap.ont.minimap2.bam >Ajap.ont.minimap2.sorted.bam
samtools index -@56 Ajap.ont.minimap2.sorted.bam
samtools view -F3844 -bSh -@ 112 -o Ajap.ont.minimap2.primaryOnly.bam Ajap.ont.minimap2.sam
samtools sort -@112 Ajap.ont.minimap2.primaryOnly.bam >Ajap.ont.minimap2.primaryOnly.sorted.bam
samtools index -@56 Ajap.ont.minimap2.primaryOnly.sorted.bam
minimap2 -ax map-ont -t 112 -o Aful.ont.minimap2.sam  Aful.genome.LG.fa Aful.ont.fastq.gz
samtools view -F4 -bSh -@ 112 -o Aful.ont.minimap2.bam Aful.ont.minimap2.sam
samtools sort -@112 Aful.ont.minimap2.bam >Aful.ont.minimap2.sorted.bam
samtools index -@56 Aful.ont.minimap2.sorted.bam
samtools view -F3844 -bSh -@ 112 -o Aful.ont.minimap2.primaryOnly.bam Aful.ont.minimap2.sam
samtools sort -@112 Aful.ont.minimap2.primaryOnly.bam >Aful.ont.minimap2.primaryOnly.sorted.bam
samtools index -@56 Aful.ont.minimap2.primaryOnly.sorted.bam
samtools view -@ 112 Ajap.ont.minimap2.sorted.bam |cut -f 1 > Ajap.ont.minimap2.AllReadName
sort Ajap.ont.minimap2.AllReadName |uniq -iu >Ajap.ont.minimap2.UniqueReadName
grep -w -f Ajap.ont.minimap2.UniqueReadName Ajap.ont.minimap2.sam >Ajap.ont.minimap2.unique.NoHeader.sam
grep "^@" Ajap.ont.minimap2.sam > Ajap.ont.minimap2.sam.header
cat Ajap.ont.minimap2.sam.header Ajap.ont.minimap2.unique.NoHeader.sam >Ajap.ont.minimap2.unique.sam
samtools view -@ 112 Aful.ont.minimap2.sorted.bam |cut -f 1 > Aful.ont.minimap2.AllReadName
sort Aful.ont.minimap2.AllReadName |uniq -iu >Aful.ont.minimap2.UniqueReadName
grep -w -f Aful.ont.minimap2.UniqueReadName Aful.ont.minimap2.sam >Aful.ont.minimap2.unique.NoHeader.sam
grep "^@" Aful.ont.minimap2.sam > Aful.ont.minimap2.sam.header
cat Aful.ont.minimap2.sam.header Aful.ont.minimap2.unique.NoHeader.sam >Aful.ont.minimap2.unique.sam

#Cross mapping
bwa mem -t 56 -M ../../../01_data/Ajap.genome.LG.fa ./Aful.illumina.clean.R1.fq.gz ./Aful.illumina.clean.R2.fq.gz -o Aful.AjapRef.illumina.sam
samtools view -bSh -@ 112 -o Aful.AjapRef.illumina.bam Aful.AjapRef.illumina.sam
samtools sort -@112 Aful.AjapRef.illumina.bam >Aful.AjapRef.illumina.sorted.bam
samtools index -@56 Aful.AjapRef.illumina.sorted.bam
samtools flagstats -@112 Aful.AjapRef.illumina.sorted.bam

bwa mem -t 56 -M ../../../01_data/Aful.genome.LG.fa ./Ajap.illumina.clean.R1.fq.gz ./Ajap.illumina.clean.R2.fq.gz -o Ajap.AfulRef.illumina.sam
samtools view -bSh -@ 112 -o Ajap.AfulRef.illumina.bam Ajap.AfulRef.illumina.sam
samtools sort -@112 Ajap.AfulRef.illumina.bam >Ajap.AfulRef.illumina.sorted.bam
samtools index -@56 Ajap.AfulRef.illumina.sorted.bam
samtools flagstats -@112 Ajap.AfulRef.illumina.sorted.bam