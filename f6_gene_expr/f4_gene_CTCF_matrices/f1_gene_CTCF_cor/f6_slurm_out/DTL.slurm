#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DTL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DTL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DTNA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DTNB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DTNBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DTWD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DTWD2 
