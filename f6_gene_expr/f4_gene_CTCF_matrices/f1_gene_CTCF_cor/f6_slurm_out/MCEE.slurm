#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MCEE.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MCEE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MCEMP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MCF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MCF2L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MCF2L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MCFD2 
