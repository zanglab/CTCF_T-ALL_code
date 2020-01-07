#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF485.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF485 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF486 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF487 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF488 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF490 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF491 
