#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PRPF38A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PRPF38A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRPF38B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRPF39 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRPF4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRPF40A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRPF40B 
