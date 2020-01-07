#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DCTD.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DCTD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCTN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCTN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCTN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCTN4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCTN5 
