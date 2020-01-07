#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DQX1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DQX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DRAM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DRAM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DRAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DRAXIN 
