#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GPR61.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR61 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR62 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR63 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR65 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR68 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR75 
