#!/bin/bash

#Peter Laurin
#script to start batch array for pairwise distances, concat segments
#

sbatch ed_calc_worker.s 
for i in {1..33}; do cat ed$i.txt >> distances_total.txt; done;
rm ed*.txt
