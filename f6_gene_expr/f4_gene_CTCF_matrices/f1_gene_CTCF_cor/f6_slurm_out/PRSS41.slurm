#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PRSS41.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PRSS41 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRSS42 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRSS45 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRSS48 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRSS50 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRSS53 
