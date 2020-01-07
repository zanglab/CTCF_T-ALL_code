#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SCARA5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SCARA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCARB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCARB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCARF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCARF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCART1 
