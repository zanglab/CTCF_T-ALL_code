#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC14A2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC14A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC15A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC15A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC15A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC15A4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC15A5 
