#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TGFBI.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TGFBI 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TGFBR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TGFBR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TGFBR3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TGFBR3L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TGFBRAP1 
