#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NUP88.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NUP88 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUP93 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUP98 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUPL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUPR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUS1 
