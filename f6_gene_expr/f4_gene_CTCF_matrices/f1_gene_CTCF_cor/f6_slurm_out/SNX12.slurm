#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SNX12.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SNX12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNX13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNX14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNX15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNX16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNX17 
