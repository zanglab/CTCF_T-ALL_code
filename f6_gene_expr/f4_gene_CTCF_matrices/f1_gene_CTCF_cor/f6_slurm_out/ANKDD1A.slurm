#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ANKDD1A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKDD1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKDD1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKEF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKFN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKFY1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKH 
