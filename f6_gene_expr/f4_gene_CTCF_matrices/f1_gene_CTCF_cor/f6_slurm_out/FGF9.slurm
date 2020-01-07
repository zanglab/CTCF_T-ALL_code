#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FGF9.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FGF9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGFBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGFBP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGFBP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGFR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGFR1OP 
