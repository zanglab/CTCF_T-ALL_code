#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o P2RY10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i P2RY10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i P2RY11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i P2RY12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i P2RY13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i P2RY14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i P2RY2 
