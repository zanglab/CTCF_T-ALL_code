#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ETV3L.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ETV3L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ETV4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ETV5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ETV6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ETV7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EVA1A 
