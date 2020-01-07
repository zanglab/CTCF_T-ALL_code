#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LINC00339.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LINC00339 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LINC00431 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LINC00475 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LINC00483 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LINC00634 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LINC00643 
