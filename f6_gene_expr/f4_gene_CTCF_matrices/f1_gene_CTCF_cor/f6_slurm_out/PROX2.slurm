#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PROX2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PROX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PROZ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRPF18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRPF19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRPF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRPF31 
