#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o POLD2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i POLD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLDIP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLDIP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLE 
