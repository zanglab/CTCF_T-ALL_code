#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CEBPE.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CEBPE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEBPG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEBPZ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEBPZOS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CECR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CECR7 
