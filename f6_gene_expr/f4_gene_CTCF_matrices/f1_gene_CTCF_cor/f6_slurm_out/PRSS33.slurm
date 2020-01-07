#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PRSS33.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PRSS33 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRSS35 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRSS36 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRSS37 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRSS38 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRSS3P2 
