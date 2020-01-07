#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CLTC.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CLTC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLTCL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLU 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLUAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLUH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLUHP3 
