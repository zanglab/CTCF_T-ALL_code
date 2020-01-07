#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC31A1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC31A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC31A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC32A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC33A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC34A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC34A2 
