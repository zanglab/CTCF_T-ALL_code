#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o POLR2K.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i POLR2K 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLR2L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLR2M 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLR3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLR3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLR3C 
