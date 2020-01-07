#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DOPEY1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DOPEY1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DOPEY2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DOT1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPAGT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPCD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPCR1 
