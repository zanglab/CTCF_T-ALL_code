#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NBPF3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NBPF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NBPF4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NBPF6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NBPF7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NBPF8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NBPF9 
