#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GP6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GP9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPA33 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPAA1 
