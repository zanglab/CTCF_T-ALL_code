#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SESN1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SESN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SESN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SESN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SESTD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SET 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SETBP1 
