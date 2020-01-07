#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GPATCH1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GPATCH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPATCH11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPATCH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPATCH2L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPATCH3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPATCH4 
