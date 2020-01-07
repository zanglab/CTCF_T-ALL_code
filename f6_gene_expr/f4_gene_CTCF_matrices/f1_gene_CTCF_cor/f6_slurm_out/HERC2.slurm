#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HERC2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HERC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HERC2P10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HERC2P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HERC2P3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HERC2P4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HERC2P7 
