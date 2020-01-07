#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CYP7B1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP7B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP8B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYR61 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYSLTR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYSLTR2 
