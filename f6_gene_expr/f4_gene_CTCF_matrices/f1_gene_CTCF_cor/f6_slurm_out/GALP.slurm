#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GALP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GALP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GALR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GALR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GALR3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GALT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GAMT 
