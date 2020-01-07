#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF85.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF85 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF850 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF852 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF853 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF860 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF862 
