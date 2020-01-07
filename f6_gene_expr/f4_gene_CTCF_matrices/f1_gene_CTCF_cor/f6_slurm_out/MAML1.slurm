#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MAML1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MAML1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAML2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAML3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAMLD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAMSTR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAN1A1 
