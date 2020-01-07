#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o UTP3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i UTP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UTP6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UTRN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UTS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UTS2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UTS2R 
