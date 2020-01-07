#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CBL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CBL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CBLB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CBLC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CBLL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CBLN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CBLN2 
