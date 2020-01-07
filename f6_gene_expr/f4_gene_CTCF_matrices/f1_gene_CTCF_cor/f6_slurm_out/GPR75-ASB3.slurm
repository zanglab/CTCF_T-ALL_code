#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GPR75-ASB3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR75-ASB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR78 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR82 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR83 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR84 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR85 
