# Metagenome-CoAssembly
### 1. Concatenate metagenomes for co-assembly 
```
cat MG_A_1.fastq MG_B_1.fastq MG_C_1.fastq MG_D_1.fastq MG_E_1.fastq > MG_1.fastq
cat MG_A_2.fastq MG_B_2.fastq MG_C_2.fastq MG_D_2.fastq MG_E_2.fastq > MG_2.fastq
```
### 2. Assemble QC'd and Trimmed metagenomic reads
```
export reads=/scratch/mpavia/MG
export output=/scratch/mpavia/assembly
megahit -1 $reads/MG_1.fastq -2 $reads/MG_2.fastq --presets meta-large -t 40 -m 0.9 -o $output/assembly
```
### 3. Quality control contigs
```
export reads=/scratch/mpavia/MG
export contigs=/scratch/mpavia/assembly
bowtie2-build -f $contigs/contigs.fa $contigs/MG.build
bowtie2 --very-sensitive-local -x $contigs/MG.build -1 $reads/MG_1.fastq -2 $reads/MG_2.fastq -S $reads/aln-MG -p 28
samtools view -@ 28 -b -S -o $reads/aln-MG.bam $reads/aln-MG
samtools sort -@ 28 $reads/aln-MG.bam -o $reads/aln-MG.bam.sorted
samtools index -@ 28 $reads/aln-MG.bam.sorted
grep ">" $contigs/contigs.fa | sed -e 's/ /\t/g' | cut -d$'\t' -f1,4 | sed 's/>//g' | sed 's/len=//g' > $contigs/contigs.lengths
genomeCoverageBed -d -split -ibam $reads/aln-MG.bam.sorted -g $contigs/contigs.fa | contig_cov_from_bed > $contigs/contig.non_zero_cov_perc.txt
awk '$3 >= 90 { print $2,$1,$3 }' $contigs/contigs.non_zero_cov_perc.txt | sort -nr | awk 'BEGIN { OFS=" " } { print $2,$1,$3 }' > $contigs/MG_select_contigs.txt
grep_ids $contigs/MG_select_contigs.txt $contigs/contigs.fa > $contigs/QC.contigs.fa
```
### 4. Determine quality of assembly
```
quast QC.contigs.fa -o quast_output
```
### 5. Build abundance files for binning
```
export reads=/scratch/mpavia/MG
export coverage=/scratch/mpavia/coverage
export contigs=/scratch/mpavia/assembly
bowtie2-build -f $contigs/QC.contigs.fa $coverage/QC.contigs.build
for i in $reads/*1.fastq;do
	echo mapping ${i:24:4}
	bowtie2 --very-sensitive-local -x $coverage/QC.contigs.build -1 $coverage/${i:24:4}_1.fastq -2 $coverage/${i:24:4}_2.fastq -S $coverage/aln-${i:24:4} -p 28
	echo sorting and indexing ${i:24:4}
	samtools view -@ 28 -b -S -o $coverage/aln-${i:24:4}.bam $coverage/aln-${i:24:4}	
        samtools sort -@ 28 $coverage/aln-${i:24:4}.bam -o $coverage/${i:24:4}.bam.sorted
	samtools index -@ 28 $coverage/${i:24:4}.bam.sorted
	echo removing ${i:24:4} intermediate files
	rm $coverage/aln-*
done
```
### 6. Binning with metabat
```
export bins=/scratch/mpavia/bins/metabat
export coverage=/scratch/mpavia/coverage
export contigs=/scratch/mpavia/assembly
jgi_summarize_bam_contig_depths --outputDepth $coverage/MG_depth.txt $coverage/MG_A.fastq.bam.sorted $coverage/MG_B.fastq.bam.sorted $coverage/MG_C.fastq.bam.sorted $coverage/MG_D.fastq.bam.sorted $coverage/MG_E.fastq.bam.sorted 
metabat -i $contigs/QC.contigs.fa -o $bins/metabat -a $coverage/MG_depth.txt -m 2000 -t 28
```
### 7. Binning with CONCOCT
```
export bins=/scratch/mpavia/bins/concoct
export coverage=/scratch/mpavia/coverage
export contigs=/scratch/mpavia/assembly
#depth files for concoct
awk '{ for (i=3;i<=NF;i+=2) $i="" }1' $coverage/MG_depth.txt > $coverage/MG_remove.abund
awk '{$2=""; print $0}' $coverage/MG_remove.abund > $coverage/MG_depth.abund
concoct --coverage_file $coverage/MG_depth.abund --composition_file $contigs/QC.contigs.fa --length_threshold 2000 -b $concoct/concoct
```
### 8. Binning with maxbin
```
export bins=/scratch/mpavia/bins/maxbin
export coverage=/scratch/mpavia/coverage
export contigs=/scratch/mpavia/assembly
#depth filess for maxbin
for a in $coverage/*.bam.sorted;do
     jgi_summarize_bam_contig_depths --outputDepth $coverage/${i:24:4}.depth.txt $coverage/${i:24:4}.bam.sorted
     cut -d$'\t' -f1,3 $coverage/${i:24:4}.depth.txt > $coverage/${i:24:4}.maxbin.abund
done
run_MaxBin.pl -contig $contigs/QC.contigs.fa -out $bins/maxbin107 -min_contig_length 2000 -thread 28 -abund $coverage/MG_A.maxbin.abund -abund1 $coverage/MG_B.maxbin.abund -abund2 $coverage/MG_C.maxbin.abund -abund3 $coverage/MG_D.maxbin.abund -abund4 $coverage/MG_E.maxbin.abund 
--markerset 107
run_MaxBin.pl -contig $contigs/QC.contigs.fa -out $bins/maxbin40 -min_contig_length 2000 -thread 28 -abund $coverage/MG_A.maxbin.abund -abund1 $coverage/MG_B.maxbin.abund -abund2 $coverage/MG_C.maxbin.abund -abund3 $coverage/MG_D.maxbin.abund -abund4 $coverage/MG_E.maxbin.abund 
--markerset 40
```
### 9. Post-processing with DASTool
```
export bins=/scratch/mpavia/bins
export contigs=/scratch/mpavia/assembly
/home/mpavia/Fasta_to_Scaffolds2Bin.sh -e fasta -i $bins/maxbin > /home/mpavia/bins/DASTool/maxbin.scaffolds2bin.tsv
/home/mpavia/Fasta_to_Scaffolds2Bin.sh -e fasta -i $bins/cococt > /home/mpavia/bins/DASTool/concoct.scaffolds2bin.tsv
/home/mpavia/Fasta_to_Scaffolds2Bin.sh -e fasta -i $bins/metabat107 > /home/mpavia/bins/DASTool/metabat107.scaffolds2bin.tsv
/home/mpavia/Fasta_to_Scaffolds2Bin.sh -e fasta -i $bins/metabat40 > /home/mpavia/bins/DASTool/metabat40.scaffolds2bin.tsv
#run in DASTool folder
DAS_Tool -i maxbin.scaffolds2bin.tsv,concoct.scaffolds2bin.tsv,metabat107.scaffolds2bin.tsv,metabat40.scaffolds2bin.tsv -c $contigs/QC.contigs.fa -o output_DASTool --score_threshold 0.2 -t 12 --search_engine blast -l maxbin,concoct,maxbin107,maxbin40
```
### 10. Extract DASTool Bins
```
export contigs=/scratch/mpavia/assembly
#run in DASTool folder
cat output_DASTool_DASTool_scaffolds2bin.txt | tr "\\t" "," > output_DASTool_DASTool_scaffolds2bin.csv
extract_fasta_bins.py $contigs/QC.contigs.fa output_DASTool_DASTool_scaffolds2bin.csv --output_path /home/mpavia/bins/DASTool
```
### 11. Check Quality of MAGs
```
checkm lineage_wf /home/mpavia/bins/DASTool /home/mpavia/bins/checkm --no_refinement -x fa --threads 12 --pplacer_threads 2
checkm qa /home/mpavia/bins/checkm/lineage.ms /home/mpavia/bins/checkm --out_format 2 --threads 1 --tab_table -f /home/mpavia/bins/checkm_output.tsv
```


