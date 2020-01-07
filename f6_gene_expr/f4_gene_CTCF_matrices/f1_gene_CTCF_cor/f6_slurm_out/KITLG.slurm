#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KITLG.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KITLG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIZ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLC2 
