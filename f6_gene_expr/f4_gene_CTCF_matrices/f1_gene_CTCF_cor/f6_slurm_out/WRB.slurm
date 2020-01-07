#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o WRB.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i WRB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WRN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WRNIP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WSB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WSB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WSCD1 
