#!/bin/bash
echo "script started"
subrange=10000
chr_id=22
for alpha in 150  # alpha
do
    for delta in 2 5 8 # delta (1% 3% 5%)
    do
        qsub -v "a=${alpha}, d=${delta}, s=${subrange}, id=${chr_id}" run_phoenix.pbs &
           
    done
done
wait
echo "All Scripts for ED haplotype aware evaluated"