#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RUFY3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RUFY3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RUFY4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RUNDC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RUNDC3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RUNDC3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RUNX1 
