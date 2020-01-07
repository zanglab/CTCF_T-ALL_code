#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TNFRSF25.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFRSF25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFRSF4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFRSF6B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFRSF8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFRSF9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFSF10 
