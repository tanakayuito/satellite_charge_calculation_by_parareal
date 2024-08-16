#!/bin/bash
#============ Slurm Options ===========
#SBATCH -p gr10451a
#SBATCH -t 1:00:00
#SBATCH --rsc p=8:t=1:c=1:m=4G
#SBATCH -o std.out
#SBATCH -e err.out

#============ Shell Script ============
set -x

srun -n 8 ./parareal_mpi_exe test_params.inp
srun -n 1 ./seq1_exe test_params.inp seq1_res.dat
srun -n 1 ./seq2_exe test_params.inp seq2_res.dat

python3.11 create_fig_seq_results.py test_params.inp seq1_res.dat seq2_res.dat
python3.11 create_fig_parareal_results.py test_params.inp seq1_res.dat

