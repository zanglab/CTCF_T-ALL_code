#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CEP68.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CEP68 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEP70 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEP72 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEP76 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEP78 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEP83 
