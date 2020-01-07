#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o WISP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i WISP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WISP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WISP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WIZ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WLS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WNK1 
