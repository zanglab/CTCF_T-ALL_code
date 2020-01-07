#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TREX2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TREX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRHDE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRHR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRIB1 
