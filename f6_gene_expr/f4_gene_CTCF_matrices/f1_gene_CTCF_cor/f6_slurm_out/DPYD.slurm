#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DPYD.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DPYD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPYS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPYSL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPYSL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPYSL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPYSL5 
