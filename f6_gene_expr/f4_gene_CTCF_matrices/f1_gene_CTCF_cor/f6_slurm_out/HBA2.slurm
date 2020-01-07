#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HBA2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HBA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HBB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HBBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HBD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HBE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HBEGF 
