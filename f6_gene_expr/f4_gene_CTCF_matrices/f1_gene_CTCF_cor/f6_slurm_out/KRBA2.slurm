#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KRBA2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KRBA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRBOX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRBOX4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRCC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KREMEN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KREMEN2 
