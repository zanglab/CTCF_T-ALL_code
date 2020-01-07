#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LY6H.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LY6H 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LY6K 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LY75 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LY75-CD302 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LY86 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LY9 
