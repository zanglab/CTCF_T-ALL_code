#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DEFA8P.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DEFA8P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DEFB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DEFB103B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DEFB108B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DEFB119 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DEFB122 
