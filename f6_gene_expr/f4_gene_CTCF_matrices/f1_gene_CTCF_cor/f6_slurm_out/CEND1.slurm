#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CEND1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CEND1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CENPA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CENPB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CENPBD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CENPBD1P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CENPC 
