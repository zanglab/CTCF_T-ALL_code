#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NBPF13P.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NBPF13P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NBPF14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NBPF15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NBPF20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NBPF22P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NBPF25P 
