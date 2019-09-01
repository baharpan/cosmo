#!/bin/bash

mkdir -p kmc_temp
ls -1 --color=no *.fastq |xargs -l -i echo "~/kmc -b -fq -ci$c -k32 -cs250 {} {}.kmc kmc_temp" >kmercount.sh
source kmercount.sh
ls -1 --color=no *.fastq |xargs -l -i echo "~/kmc_tools sort {}.kmc {}.kmc.sorted " >kmercountsort.sh
source kmercountsort.sh
ls -1 --color=no *.fastq |xargs -l -i echo "{}.kmc.sorted" > filtered_kmc2_list

file_lines=($(wc -l test.fastq))
n_seqs="$((${file_lines}/4))"

numactl --interleave=all  ../cosmo-pack -k filtered_kmc2_list
/usr/bin/time -v ../bubbles_matrix -n $n_seqs -k 32 -i test.fastq
/usr/bin/time -v ../bubbles -c 3 -u 2 -b 1000  -k 32 > output
