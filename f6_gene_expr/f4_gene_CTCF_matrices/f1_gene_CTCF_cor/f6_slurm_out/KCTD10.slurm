#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KCTD10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KCTD10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCTD11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCTD12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCTD13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCTD14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCTD15 
