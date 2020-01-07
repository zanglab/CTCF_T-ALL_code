#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PCDHGA8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHGA8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHGA9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHGB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHGB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHGB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHGB4 
