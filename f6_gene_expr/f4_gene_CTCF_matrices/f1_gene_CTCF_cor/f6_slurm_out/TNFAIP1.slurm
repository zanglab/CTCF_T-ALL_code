#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TNFAIP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFAIP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFAIP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFAIP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFAIP6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFAIP8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFAIP8L1 
