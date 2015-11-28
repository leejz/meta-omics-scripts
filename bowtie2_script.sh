#!/bin/bash -l
#PBS -q regular
#PBS -l mppwidth=120
#PBS -l walltime=23:00:00
#PBS -N bowtie_2_MT_k29
#PBS -o $PBS_JOBNAME-$PBS_JOBID.out
#PBS -j oe
#PBS -V

cd $PBS_O_WORKDIR
module load bowtie2
PAIREDFILES=./MT_source/*.1.fastq 

for i in $PAIREDFILES
do
  aprun -b -n 1 bowtie2 -x ./k29.scaffolds.bowtie2.index -1 ./$i -2 ./${i%.1.fastq}.2.fastq -S ./${i%.1.fastq}.k29.scaffolds.mapping.sam &
done
wait
