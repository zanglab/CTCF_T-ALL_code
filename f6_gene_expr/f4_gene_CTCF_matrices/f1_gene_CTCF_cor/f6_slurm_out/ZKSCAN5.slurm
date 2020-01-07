#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZKSCAN5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZKSCAN5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZKSCAN7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZKSCAN8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZMAT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZMAT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZMAT3 
