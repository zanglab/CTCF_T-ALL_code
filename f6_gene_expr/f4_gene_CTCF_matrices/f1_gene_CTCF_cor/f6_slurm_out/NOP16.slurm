#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NOP16.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NOP16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOP56 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOP58 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOP9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOS1 
