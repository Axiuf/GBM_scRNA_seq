#!/bin/bash

for win_len in 11 21 31 51 81 101
do
	Rscript ./src/7_GBM_total_merge_InferCNV_test.R $win_len &
done

wait
