#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KCNT2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNU1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNV1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCNV2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCTD1 
