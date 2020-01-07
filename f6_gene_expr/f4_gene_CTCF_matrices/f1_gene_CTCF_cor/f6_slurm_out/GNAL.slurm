#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GNAL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GNAL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNAO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNAQ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNAS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNAT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GNAT2 
