#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMEFF1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEFF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEFF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM100 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM101 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM102 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMEM104 
