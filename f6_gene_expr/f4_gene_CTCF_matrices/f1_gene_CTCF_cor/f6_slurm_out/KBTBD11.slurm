#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KBTBD11.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KBTBD11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KBTBD11-OT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KBTBD12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KBTBD13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KBTBD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KBTBD3 
