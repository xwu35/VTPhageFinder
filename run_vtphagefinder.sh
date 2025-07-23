#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH -J run_VTPhageFinder
#SBATCH --output=slurm-%x.%j.out
#SBATCH --error=slurm-%x.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xiaofen

eval "$(conda shell.bash hook)"
conda activate snakemake

python VTPhageFinder.py \
--reads_dir test_data/sequences \
--sample_info test_data/sample_info.txt \
--output_dir FpVT_output \
--reference_genome resources/Fp22_genome/fp22_assembly.fasta \
--prophage_region resources/Fp22_genome/fp22_prophage_region.bed
