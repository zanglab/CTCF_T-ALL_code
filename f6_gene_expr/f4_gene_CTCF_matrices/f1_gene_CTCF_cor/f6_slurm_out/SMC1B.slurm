#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SMC1B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SMC1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMC5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SMC6 
