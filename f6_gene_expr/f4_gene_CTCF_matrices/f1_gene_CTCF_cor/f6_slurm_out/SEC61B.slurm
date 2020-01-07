#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SEC61B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SEC61B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEC61G 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEC62 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEC63 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SECISBP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SECISBP2L 
