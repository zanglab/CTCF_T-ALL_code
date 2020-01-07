#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DSTN.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DSTN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DSTNP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DSTYK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DTD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DTD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DTHD1 
