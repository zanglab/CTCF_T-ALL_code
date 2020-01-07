#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GPR87.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR87 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR88 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR89A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPR89B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPRASP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPRASP2 
