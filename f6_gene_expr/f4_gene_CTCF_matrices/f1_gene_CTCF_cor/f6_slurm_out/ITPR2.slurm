#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ITPR2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ITPR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITPR3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITPRIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITPRIPL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITPRIPL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ITSN1 
