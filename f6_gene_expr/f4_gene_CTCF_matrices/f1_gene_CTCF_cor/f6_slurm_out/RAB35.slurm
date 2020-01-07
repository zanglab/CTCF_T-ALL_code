#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RAB35.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB35 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB36 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB37 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB38 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB39A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB39B 
