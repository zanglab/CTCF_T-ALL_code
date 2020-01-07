#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DSCAML1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DSCAML1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DSCC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DSE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DSEL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DSG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DSG2 
