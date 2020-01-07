#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PROKR1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PROKR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PROKR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PROM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PROM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PROP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRORSD1P 
