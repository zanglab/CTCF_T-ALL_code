#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HIST1H4I.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H4I 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H4J 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H4K 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H4L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST2H2AA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST2H2AA4 
