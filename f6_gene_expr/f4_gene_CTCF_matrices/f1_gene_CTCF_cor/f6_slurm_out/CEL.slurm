#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CEL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CEL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CELA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CELA2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CELA2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CELA3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CELA3B 
