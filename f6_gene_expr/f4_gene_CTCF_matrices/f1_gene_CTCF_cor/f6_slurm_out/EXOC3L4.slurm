#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EXOC3L4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EXOC3L4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EXOC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EXOC5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EXOC6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EXOC6B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EXOC7 
