#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PCNA.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PCNA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCNAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCNP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCNT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCOLCE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCOLCE2 
