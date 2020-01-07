#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC10A1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC10A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC10A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC10A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC10A4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC10A5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC10A6 
