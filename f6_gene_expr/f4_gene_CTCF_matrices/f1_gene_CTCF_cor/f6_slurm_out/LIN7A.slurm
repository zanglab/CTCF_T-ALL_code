#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LIN7A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LIN7A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LIN7B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LIN7C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LIN9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LINC00032 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LINC00301 
