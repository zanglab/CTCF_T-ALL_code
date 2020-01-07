#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SYN2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SYN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYNC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYNCRIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYNDIG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYNDIG1L 
