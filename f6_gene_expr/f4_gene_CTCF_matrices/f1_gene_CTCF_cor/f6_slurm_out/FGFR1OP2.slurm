#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FGFR1OP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FGFR1OP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGFR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGFR3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGFR4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGFRL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGG 
