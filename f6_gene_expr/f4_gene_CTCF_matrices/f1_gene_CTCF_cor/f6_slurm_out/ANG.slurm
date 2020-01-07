#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ANG.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ANG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANGEL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANGEL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANGPT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANGPT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANGPT4 
