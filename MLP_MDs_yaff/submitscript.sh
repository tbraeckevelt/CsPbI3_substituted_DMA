#!/bin/sh

#PBS -o output_main.txt
#PBS -e error_main.txt
#PBS -l nodes=1:ppn=16
#PBS -A 2024_026
#PBS -l walltime=71:59:00

cd ${PBS_O_WORKDIR}

ml --force purge
ml psiflow/3.0.0_own

cp utils_debug.py /data/gent/vo/000/gvo00003/shared/tombraeckevelt/clusters/gpu_rome_a100/software/psiflow-3.0.0_own/psiflow_env/lib/python3.9/site-packages/.

python main.py hortense.yaml
