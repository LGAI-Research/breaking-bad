#!/bin/bash
#SBATCH --job-name=lg-scd-exp
#SBATCH --cpus-per-task=8                 
#SBATCH --mem=32G                         
#SBATCH --time=2-00:00:00
#SBATCH --nodelist=colossal 

# 1. 파일들이 있는 실제 폴더로 기준 경로 변경
#SBATCH --chdir=/home/taehui/LG_SCD/breaking-bad/paper_experiments

# (로그 파일은 원하는 위치 그대로)
#SBATCH --output=/home/taehui/LG_SCD/slurm-logs/%j.out
#SBATCH --error=/home/taehui/LG_SCD/slurm-logs/%j.err
 
mkdir -p /home/taehui/LG_SCD/slurm-logs
 
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK
 
source /home/taehui/home/taehui/miniforge3/etc/profile.d/conda.sh
conda activate tei_breaking_bad
 
# 2. --chdir로 위치를 맞췄으니 쉘 스크립트만 바로 실행
bash run_all_experiments.sh
#python run_boss_tetrad.py