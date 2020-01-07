#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TTK.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TTK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTLL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTLL10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTLL11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTLL12 
