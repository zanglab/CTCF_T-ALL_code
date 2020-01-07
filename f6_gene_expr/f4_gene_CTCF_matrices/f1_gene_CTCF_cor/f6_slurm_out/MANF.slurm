#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MANF.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MANF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MANSC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MANSC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAOA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAOB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAP10 
