#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o POP7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i POP7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POPDC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POPDC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PORCN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POSTN 
