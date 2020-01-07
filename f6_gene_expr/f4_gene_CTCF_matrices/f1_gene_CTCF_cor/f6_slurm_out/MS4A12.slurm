#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MS4A12.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MS4A12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MS4A13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MS4A14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MS4A15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MS4A18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MS4A2 
