#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TTC28.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TTC28 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTC29 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTC30A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTC30B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TTC31 
