#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RAX2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RAX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RB1CC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBAK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBAK-RBAKDN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBBP4 
