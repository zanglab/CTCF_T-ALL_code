#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HOXC6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HOXC6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HOXC8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HOXC9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HOXD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HOXD10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HOXD11 
