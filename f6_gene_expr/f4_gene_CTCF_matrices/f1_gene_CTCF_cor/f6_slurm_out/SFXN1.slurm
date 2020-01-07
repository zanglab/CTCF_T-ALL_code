#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SFXN1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SFXN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SFXN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SFXN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SFXN4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SFXN5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SGCA 
