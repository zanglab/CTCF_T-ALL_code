#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MECR.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MECR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MED1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MED10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MED11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MED12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MED12L 
