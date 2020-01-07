#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CISD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CISD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CISD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CISD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CISH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CIT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CITED1 
