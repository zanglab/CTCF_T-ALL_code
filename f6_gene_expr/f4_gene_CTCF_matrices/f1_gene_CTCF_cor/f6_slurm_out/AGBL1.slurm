#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AGBL1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AGBL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGBL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGBL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGBL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGBL5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGER 
