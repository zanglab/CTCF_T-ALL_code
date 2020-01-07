#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TEKT4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TEKT4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEKT4P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEKT5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TELO2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TENM1 
