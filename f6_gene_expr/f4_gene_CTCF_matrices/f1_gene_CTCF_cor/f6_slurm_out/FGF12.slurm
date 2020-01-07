#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FGF12.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FGF12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGF13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGF14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGF16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGF17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FGF18 
