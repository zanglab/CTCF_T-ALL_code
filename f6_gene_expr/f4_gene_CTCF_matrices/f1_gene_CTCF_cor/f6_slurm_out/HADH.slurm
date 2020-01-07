#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HADH.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HADH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HADHA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HADHB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HAGH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HAGHL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HAL 
