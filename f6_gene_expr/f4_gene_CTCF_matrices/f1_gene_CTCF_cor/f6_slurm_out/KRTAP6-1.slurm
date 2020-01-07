#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KRTAP6-1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP6-1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP7-1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP8-1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP9-4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTAP9-6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRTCAP2 
