#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PCDHGA2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHGA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHGA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHGA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHGA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHGA6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDHGA7 
