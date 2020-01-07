#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o STATH.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i STATH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STAU1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STAU2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STBD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STC2 
