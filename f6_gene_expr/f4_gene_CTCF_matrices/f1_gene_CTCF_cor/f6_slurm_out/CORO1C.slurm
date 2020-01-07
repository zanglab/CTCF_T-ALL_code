#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CORO1C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CORO1C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CORO2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CORO2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CORO6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CORO7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CORO7-PAM16 
