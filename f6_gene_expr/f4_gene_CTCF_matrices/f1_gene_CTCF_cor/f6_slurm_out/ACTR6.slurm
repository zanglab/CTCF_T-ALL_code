#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ACTR6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ACTR6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACTR8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACTRT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACTRT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACTRT3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACVR1 
