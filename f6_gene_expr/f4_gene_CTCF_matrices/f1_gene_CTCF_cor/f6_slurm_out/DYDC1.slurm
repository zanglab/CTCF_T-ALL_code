#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DYDC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DYDC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DYDC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DYM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DYNAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DYNC1H1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DYNC1I1 
