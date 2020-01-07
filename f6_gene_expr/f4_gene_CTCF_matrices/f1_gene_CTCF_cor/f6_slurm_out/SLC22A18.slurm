#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC22A18.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC22A18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC22A18AS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC22A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC22A23 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC22A24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC22A25 
