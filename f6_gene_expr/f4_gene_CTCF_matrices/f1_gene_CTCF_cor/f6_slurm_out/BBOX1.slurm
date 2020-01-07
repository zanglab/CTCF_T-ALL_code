#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BBOX1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BBOX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BBS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BBS10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BBS12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BBS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BBS4 
