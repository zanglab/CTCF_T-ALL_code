#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KLK1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KLK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLK10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLK11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLK12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLK13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLK14 
