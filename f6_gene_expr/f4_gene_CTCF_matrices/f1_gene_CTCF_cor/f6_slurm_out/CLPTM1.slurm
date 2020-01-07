#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CLPTM1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CLPTM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLPTM1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLPX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLRN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLRN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLRN3 
