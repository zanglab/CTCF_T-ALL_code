#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF442.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF442 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF443 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF444 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF445 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF446 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF449 
