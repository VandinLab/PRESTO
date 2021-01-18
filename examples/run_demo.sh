#!/bin/bash
fold_ex="exact/"
fold_par="parallel-sampling/"
fold_motifs="motifs/"
b=5
r=100
c=1.25
curr_dir=${PWD}
fold_data="test-graphs/"
delta=86400

for graph in ${fold_data}*
do
	f_name="$(basename -- $graph)"
	out_file="${f_name}_demo.txt"
	rm ${out_file}
	echo "################################" >> ${out_file} 
	echo $f_name >> ${out_file} 
	echo "################################" >> ${out_file} 
	for motif in ${fold_motifs}*
	do
		m_name="$(basename -- $motif)"
		echo "++++++++++++++++++++++++++++++++" >> ${out_file} 
		echo $m_name >> ${out_file} 
		echo "++++++++++++++++++++++++++++++++" >> ${out_file} 
		cd ${fold_ex}
		./exactBT -i:../${graph} -m:../${motif} -delta:${delta} >> ../${out_file} 2>&1
		echo "--------------------------------" >> ../${out_file} 
		cd ../${fold_par}
		./motifsamplermain -i:../${graph} -m:../${motif} -delta:${delta} -c:${c} -r:${r} -b:${b} -type:1 >> ../${out_file} 2>&1
		cd ${curr_dir}
	done
done
#done
