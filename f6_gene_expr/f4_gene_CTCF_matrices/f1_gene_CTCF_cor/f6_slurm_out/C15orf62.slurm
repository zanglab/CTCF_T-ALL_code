#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C15orf62.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C15orf62 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C15orf65 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C16orf45 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C16orf46 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C16orf54 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C16orf58 
