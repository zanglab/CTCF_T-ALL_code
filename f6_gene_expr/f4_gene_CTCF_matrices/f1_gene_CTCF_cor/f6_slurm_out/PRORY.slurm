#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PRORY.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PRORY 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PROS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PROSER1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PROSER2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PROSER3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PROX1 
