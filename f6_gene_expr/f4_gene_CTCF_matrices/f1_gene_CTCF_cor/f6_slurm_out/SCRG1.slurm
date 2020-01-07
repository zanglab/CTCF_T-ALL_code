#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SCRG1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SCRG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCRIB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCRN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCRN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCRN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCRT1 
