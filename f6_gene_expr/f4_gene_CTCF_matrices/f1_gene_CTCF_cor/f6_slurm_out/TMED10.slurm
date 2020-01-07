#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMED10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMED10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMED10P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMED11P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMED2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMED3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMED4 
