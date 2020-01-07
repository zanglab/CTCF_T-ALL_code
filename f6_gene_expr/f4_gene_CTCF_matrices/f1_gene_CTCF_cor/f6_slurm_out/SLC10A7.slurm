#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC10A7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC10A7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC11A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC11A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC12A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC12A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC12A3 
