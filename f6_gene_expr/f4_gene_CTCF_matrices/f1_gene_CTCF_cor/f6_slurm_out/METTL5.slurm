#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o METTL5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i METTL5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i METTL6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i METTL7A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i METTL7B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i METTL8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i METTL9 
