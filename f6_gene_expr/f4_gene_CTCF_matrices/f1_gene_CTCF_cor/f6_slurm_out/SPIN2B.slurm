#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SPIN2B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SPIN2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPIN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPIN4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPINK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPINK13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPINK14 
