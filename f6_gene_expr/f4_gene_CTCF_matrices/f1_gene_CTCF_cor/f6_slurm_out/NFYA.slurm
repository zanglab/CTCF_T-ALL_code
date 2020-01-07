#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NFYA.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NFYA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NFYB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NFYC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NGB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NGDN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NGEF 
