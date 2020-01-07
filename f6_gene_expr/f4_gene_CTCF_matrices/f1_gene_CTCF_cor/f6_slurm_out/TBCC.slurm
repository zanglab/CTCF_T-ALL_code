#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TBCC.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TBCC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBCCD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBCD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBCE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBCEL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBCK 
