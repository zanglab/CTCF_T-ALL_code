#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CYP19A1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP19A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP1A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP1A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP1B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP20A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYP21A1P 
