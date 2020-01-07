#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RUNX1T1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RUNX1T1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RUNX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RUNX3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RUSC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RUSC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RUVBL1 
