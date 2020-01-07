#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GAS1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GAS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GAS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GAS2L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GAS2L1P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GAS2L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GAS2L3 
