#!/bin/bash
#============ Slurm Options ===========
#SBATCH -p gr10451a
#SBATCH -t 1:00:00
#SBATCH --rsc p=4:t=1:c=1:m=16G
#SBATCH -o std.out
#SBATCH -e err.out

#============ Shell Script ============
set -x

srun -n 4 ./parareal_mpi_exe params.inp
srun -n 1 ./seq1_exe params.inp seq1_res.dat
srun -n 1 ./seq2_exe params.inp seq2_res.dat
#srun -n6 -c1 -l --multi-prog multi.conf

python3.11 create_fig_seq_results.py params.inp seq1_res.dat seq2_res.dat
python3.11 create_fig_parareal_results.py params.inp seq1_res.dat

