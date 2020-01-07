#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DOC2B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DOC2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DOC2GP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DOCK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DOCK10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DOCK11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DOCK2 
