#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NUDT8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NUDT8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUDT9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUDT9P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUFIP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUFIP2 
