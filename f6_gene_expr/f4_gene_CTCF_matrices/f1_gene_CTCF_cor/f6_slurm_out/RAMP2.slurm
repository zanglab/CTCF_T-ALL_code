#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RAMP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RAMP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAMP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RANBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RANBP10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RANBP17 
