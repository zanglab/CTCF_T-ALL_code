#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GNL3L.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GNL3L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNLY 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNMT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNPAT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNPDA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNPDA2 
