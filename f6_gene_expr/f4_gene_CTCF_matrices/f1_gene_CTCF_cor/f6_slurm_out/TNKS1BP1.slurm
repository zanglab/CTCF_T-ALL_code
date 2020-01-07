#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TNKS1BP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TNKS1BP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNKS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNMD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNNC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNNC2 
