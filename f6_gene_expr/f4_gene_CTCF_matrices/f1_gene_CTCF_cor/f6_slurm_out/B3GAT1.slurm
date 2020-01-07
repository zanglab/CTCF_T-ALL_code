#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o B3GAT1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i B3GAT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i B3GAT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i B3GAT3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i B3GLCT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i B3GNT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i B3GNT3 
