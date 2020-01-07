#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o INTS7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i INTS7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INTS8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INTS9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INTU 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INVS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IP6K1 
