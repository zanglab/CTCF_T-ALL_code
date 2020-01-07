#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF770.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF770 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF771 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF772 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF773 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF774 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF775 
