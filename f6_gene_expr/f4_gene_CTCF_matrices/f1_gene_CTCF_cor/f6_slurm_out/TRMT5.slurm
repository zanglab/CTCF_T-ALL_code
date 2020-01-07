#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TRMT5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TRMT5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRMT6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRMT61A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRMT61B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRMU 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRNAU1AP 
