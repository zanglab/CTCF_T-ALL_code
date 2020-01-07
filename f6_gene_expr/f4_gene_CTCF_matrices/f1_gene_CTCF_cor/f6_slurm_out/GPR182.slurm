#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GPR182.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR182 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR183 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR22 
