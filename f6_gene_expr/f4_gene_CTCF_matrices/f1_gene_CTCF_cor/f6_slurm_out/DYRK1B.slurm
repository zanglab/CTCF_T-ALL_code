#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DYRK1B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DYRK1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DYRK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DYRK3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DYRK4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DYSF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DYTN 
