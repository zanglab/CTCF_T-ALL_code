#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ILDR1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ILDR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ILDR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ILF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ILF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ILK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ILKAP 
