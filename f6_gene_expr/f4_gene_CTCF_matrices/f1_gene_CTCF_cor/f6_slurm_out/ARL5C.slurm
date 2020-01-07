#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ARL5C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ARL5C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARL6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARL6IP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARL6IP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARL6IP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARL6IP6 
