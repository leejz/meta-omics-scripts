#!/bin/bash -l 

echo "Running Samtools, convert sam to bam and sort"
module load samtools

SAMFILES=./*.sam
for i in $SAMFILES
do
  samtools view -bS $i > ${i%.sam}.bam &
done
wait

BAMFILES=./*.bam
for j in $BAMFILES
do
  samtools sort $j ${j%.bam}.sorted &
done
wait

echo "Running bedtools on sorted bam alignments with gff"
SORTEDFILES=./*sorted.bam
for k in $SORTEDFILES
do
  outname=$(basename $k)
  bedtools genomecov -d -ibam $k > k29.Scaffolds.${outname%.sorted.bam}.coverage.txt &
done
wait

echo "Running coverage summary script"
COVFILES=./*coverage.txt
for i in $COVFILES
do
  outname=$(basename $i)
  contig_coverage_from_perbase_gff.py -i $i -t -o ${outname%.coverage.txt}.Counts.txt & 
done
wait
