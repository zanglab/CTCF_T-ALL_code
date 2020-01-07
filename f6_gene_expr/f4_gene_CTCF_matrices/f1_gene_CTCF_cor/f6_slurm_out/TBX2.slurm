#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TBX2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TBX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBX20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBX21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBX22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBX3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBX4 
