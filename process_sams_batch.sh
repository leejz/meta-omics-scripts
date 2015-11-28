#!/bin/bash -l
#PBS -q regular
#PBS -l mppwidth=192
#PBS -l walltime=11:00:00
#PBS -N k29_MT_process_sam
#PBS -o $PBS_JOBNAME-$PBS_JOBID.out
#PBS -j oe
#PBS -V


#echo "Running Samtools, convert sam to bam and sort"
cd $PBS_O_WORKDIR
module load samtools
echo "Generating bam..."
echo "Started `date`"
SAMFILES=./*.sam
for i in $SAMFILES
do
  aprun -b -n 1 -q samtools view -bS $i > ${i%.sam}.bam &
done
wait

echo "Sorting..."
echo "Started `date`"
BAMFILES=./*.bam
for j in $BAMFILES
do
  aprun -b -n 1 samtools sort $j ${j%.bam}.sorted &
done
wait

echo "computing bed"
echo "Started `date`"

echo "Running bedtools on sorted bam alignments with gff"
SORTEDFILES=./*sorted.bam.bam
for k in $SORTEDFILES
do
  outname=$(basename $k)
  aprun -b -n 1 -q bedtools coverage -d -abam $k -b k29.MT.fa.gff > k29.MT.${outname%.sorted.bam.bam}.coverage.txt &
done
wait

echo "Running coverage summary script"
echo "Started `date`"
COVFILES=./*coverage.txt
for i in $COVFILES
do
  outname=$(basename $i)
  aprun -b -n 1 contig_coverage_from_perbase_gff.py -i $i -t -o ${outname%.coverage.txt}.Counts.txt & 
done
wait

echo "Done! at `date`"
