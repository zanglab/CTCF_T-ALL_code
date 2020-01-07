#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HSPA14.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HSPA14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSPA1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSPA1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSPA1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSPA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSPA4 
