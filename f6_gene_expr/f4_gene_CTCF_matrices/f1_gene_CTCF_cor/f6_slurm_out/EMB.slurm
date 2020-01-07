#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EMB.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EMB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EMBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EMC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EMC10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EMC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EMC3 
