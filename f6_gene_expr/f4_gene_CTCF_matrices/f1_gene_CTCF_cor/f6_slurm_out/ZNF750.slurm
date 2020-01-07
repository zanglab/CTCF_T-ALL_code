#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF750.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF750 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF75A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF75D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF76 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF761 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF763 
