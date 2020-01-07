#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TCERG1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TCERG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCERG1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCF12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCF15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCF19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCF20 
