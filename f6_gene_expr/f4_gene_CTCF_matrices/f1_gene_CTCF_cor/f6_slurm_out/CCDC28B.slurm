#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CCDC28B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC28B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC30 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC33 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC34 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC36 
