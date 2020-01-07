#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HSPA4L.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HSPA4L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSPA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSPA6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSPA7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSPA8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSPA9 
