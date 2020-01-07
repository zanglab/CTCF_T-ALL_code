#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PCDHB8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHB8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHB9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHGA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHGA10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHGA11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHGA12 
