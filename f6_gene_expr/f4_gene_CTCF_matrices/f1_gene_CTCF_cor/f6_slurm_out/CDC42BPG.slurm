#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CDC42BPG.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CDC42BPG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDC42EP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDC42EP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDC42EP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDC42EP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDC42EP5 
