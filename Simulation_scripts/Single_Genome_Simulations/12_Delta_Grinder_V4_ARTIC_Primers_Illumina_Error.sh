#!/bin/bash
#SBATCH --job-name=Grinder_sim_V4_ARTIC_Illumina_reads_Error
#SBATCH --partition=glenn_p
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=20gb
#SBATCH --export=NONE
#SBATCH --time=90:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=email@uga.edu
#SBATCH --mail-type=BEGIN,END,FAIL

module load Grinder/0.5.4-GCCcore-8.3.0-Perl-5.30.0

cd /scratch/mandyh/WISER_InSilico/Amplicon

#define directory
genomes='/scratch/mandyh/WISER_InSilico/genomes'
input='/scratch/mandyh/WISER_InSilico/Amplicon/V4_ARTIC_Primers'
output='/scratch/mandyh/WISER_InSilico/Amplicon/sim_V4_Delta_Illumina_reads_error'
final_output='/scratch/mandyh/WISER_InSilico/Amplicon/final_Illumina_error_reads'

for i in {1..99}
do
for x in {150,250,300}
do
for y in {1,2,3,4,5,6,7,8,9,10,15,20,25,30,35,40,45,50,75,100}
do
for z in {1..5}
do
echo "i: $i, x: $x, y: $y, z: $z"
grinder -reference_file $genomes/Delta_gisaid_hcov-19_GA_2021_06_30_14.fasta -forward_reverse $input/V4_POS$i\.fas -coverage_fold $y -read_dist $x -insert_dist 350 -mate_orientation FR -mutation_dist linear 0.1 0.2 -qual_levels 30 10 -fastq_output 1 -output_dir $output/ -base_name Delta_POS$i\_PE$x\_$y\X_IE_$z
done
done
done
done

for x in {150,250,300}
do
for y in {1,2,3,4,5,6,7,8,9,10,15,20,25,30,35,40,45,50,75,100}
do
for z in {1..5}
do
echo "x: $x, y: $y, z: $z"
cat $output/Delta_POS*_PE$x\_$y\X_IE_$z\-reads.fastq > $final_output/Delta_V4_PE$x\_$y\X_IE_$z\-reads.fastq
done
done
done