#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GNA13.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GNA13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNA14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNA15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNAI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNAI2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNAI3 
