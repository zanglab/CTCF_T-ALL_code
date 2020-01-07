#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SERPINA2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPINA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPINA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPINA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPINA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPINA6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPINA7 
