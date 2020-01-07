#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NPAS2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NPAS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPAS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPAS4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPAT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPBWR1 
