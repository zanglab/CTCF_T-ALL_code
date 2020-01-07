#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC40A1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC40A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC41A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC41A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC41A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC43A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC43A2 
