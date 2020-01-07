#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SCN9A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SCN9A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCNM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCNN1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCNN1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCNN1D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCNN1G 
