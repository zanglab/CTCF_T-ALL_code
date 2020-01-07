#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PPIF.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PPIF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPIG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPIH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPIL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPIL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPIL3 
