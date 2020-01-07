#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PM20D1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PM20D1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PM20D2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PMAIP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PMCH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PMCHL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PMCHL2 
