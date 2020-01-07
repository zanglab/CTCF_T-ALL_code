#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FANCF.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FANCF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FANCG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FANCI 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FANCL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FANCM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FANK1 
