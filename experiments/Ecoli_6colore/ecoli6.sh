set -xe 
mkdir -p kmc_temp
ls -1 --color=no *.fna |xargs -l -i echo "cat {} |python3 /s/fir/c/nobackup/baharpan/git/cosmo/experiments/join_fasta_lines_stream.py {}.flat" >flatten.sh
source flatten.sh 
ls -1 --color=no *.flat |xargs -l -i echo "/s/fir/c/nobackup/baharpan/git/KMC/bin/kmc -ci0 -fa -k31 -cs300 {} {}.kmc kmc_temp" >kmercount.sh
source kmercount.sh
ls -1 --color=no *.flat |xargs -l -i echo "/s/fir/c/nobackup/baharpan/git/KMC/bin/kmc_tools sort {}.kmc {}.kmc.sorted " >kmercountsort.sh
source kmercountsort.sh
ls -1 --color=no *.flat |xargs -l -i echo "{}.kmc" > ecoli6_kmc2_list

numactl --interleave=all /bin/time -v /s/fir/c/nobackup/baharpan/git/cosmo/cosmo-pack -k ecoli6_kmc2_list
numactl --interleave=all /bin/time -v /s/fir/c/nobackup/baharpan/git/cosmo/pack-color ecoli6_kmc2_list.colors 6
/bin/time -v  /s/fir/c/nobackup/baharpan/git/cosmo/cosmo-color ecoli6_kmc2_list.packed ecoli6_kmc2_list.colors.rrr
