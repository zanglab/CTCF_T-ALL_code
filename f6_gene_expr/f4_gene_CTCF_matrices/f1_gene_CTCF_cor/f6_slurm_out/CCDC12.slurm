#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CCDC12.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC120 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC121 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC122 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC124 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCDC125 
