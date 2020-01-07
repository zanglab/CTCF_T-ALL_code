#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LRP10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LRP10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRP11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRP12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRP1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRP2BP 
