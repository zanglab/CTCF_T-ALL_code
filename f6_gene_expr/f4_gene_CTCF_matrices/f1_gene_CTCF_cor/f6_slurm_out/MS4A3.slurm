#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MS4A3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MS4A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MS4A4A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MS4A5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MS4A6A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MS4A6E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MS4A7 
