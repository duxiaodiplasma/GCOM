#!/bin/bash
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=duxiaodi@fusion.gat.com
#SBATCH --job-name=const_mu
#SBATCH --output=/home/duxiaodi/GCOM/GCOM_v2/RESONANCE/JOB/30_const_mu.job
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=10000

ipython /home/duxiaodi/GCOM/GCOM_v2/RESONANCE/const_mu.py
