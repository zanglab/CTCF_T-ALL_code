#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ROBO4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ROBO4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ROCK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ROCK1P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ROCK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ROGDI 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ROM1 
