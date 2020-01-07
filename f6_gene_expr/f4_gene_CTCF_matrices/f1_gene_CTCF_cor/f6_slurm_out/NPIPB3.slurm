#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NPIPB3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NPIPB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPIPB4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPIPB5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPIPB6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPIPB8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPIPB9 
