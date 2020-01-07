#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OXGR1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OXGR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OXLD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OXNAD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OXR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OXSM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OXSR1 
