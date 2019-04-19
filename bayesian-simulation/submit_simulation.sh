#!/bin/sh
#
#
cd ~/Linbo/odds/final_final_final_simulation
seed=14

for sample_size in 50 100 200 400 800 1600
do
	qsub ./call_simulation.sh $sample_size $seed
done

