#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ARID3B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ARID3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARID3C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARID4A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARID4B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARID5A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARID5B 
