#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NPIPA3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NPIPA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPIPA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPIPA7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPIPA8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPIPB11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPIPB15 
