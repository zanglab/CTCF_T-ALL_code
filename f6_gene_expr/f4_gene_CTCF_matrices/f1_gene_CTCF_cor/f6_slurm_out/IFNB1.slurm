#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IFNB1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IFNB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFNE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFNG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFNGR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFNGR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IFNK 
