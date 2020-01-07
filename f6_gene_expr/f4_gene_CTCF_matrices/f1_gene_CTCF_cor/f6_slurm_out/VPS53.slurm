#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o VPS53.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i VPS53 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VPS54 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VPS72 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VPS8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VPS9D1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VRK1 
