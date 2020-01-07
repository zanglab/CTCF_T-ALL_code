#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NBL1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NBL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NBN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NBPF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NBPF10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NBPF11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NBPF12 
