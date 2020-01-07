#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF467.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF467 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF468 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF469 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF470 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF471 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF473 
