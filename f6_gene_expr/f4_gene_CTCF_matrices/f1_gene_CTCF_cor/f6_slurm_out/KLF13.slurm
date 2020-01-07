#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KLF13.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KLF13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLF14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLF15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLF16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLF17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLF2 
