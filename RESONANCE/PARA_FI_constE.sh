#!/bin/bash 
#SBATCH --mail-type=END,FAIL 
#SBATCH --mail-user=duxiaodi@fusion.gat.com 
#SBATCH --job-name=PARA_FI.sh
#SBATCH --output=/home/duxiaodi/GCOM/GCOM_v2/RESONANCE/OUT/FI.job
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12   
#SBATCH --time=23:59:59
#SBATCH --mem-per-cpu=10000 

ipython /home/duxiaodi/GCOM/GCOM_v2/RESONANCE/PARA_FI_constE.py
