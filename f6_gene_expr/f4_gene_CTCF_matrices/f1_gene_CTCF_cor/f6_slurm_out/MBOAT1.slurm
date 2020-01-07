#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MBOAT1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MBOAT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MBOAT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MBOAT4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MBOAT7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MBP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MBTD1 
