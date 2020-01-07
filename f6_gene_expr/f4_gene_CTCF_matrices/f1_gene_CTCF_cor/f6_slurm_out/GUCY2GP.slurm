#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GUCY2GP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GUCY2GP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GUF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GUK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GULP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GUSB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GUSBP1 
