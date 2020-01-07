#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LINC00672.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LINC00672 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LINC00674 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LINC00854 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LINC00869 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LINC00888 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LINC00933 
