#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ACSS3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ACSS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACTA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACTA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACTB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACTBL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACTC1 
