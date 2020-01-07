#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PRAMEF17.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PRAMEF17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRAMEF18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRAMEF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRAMEF20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRAMEF5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRAMEF6 
