#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GIPC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GIPC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GIPC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GIPC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GIPR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GIT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GIT2 
