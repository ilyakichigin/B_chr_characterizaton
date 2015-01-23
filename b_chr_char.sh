#!/bin/bash
IFS=$'\r\n' GLOBIGNORE='*' :; fastq_array=($(ls *.fastq | cat))
IFS=$'\r\n' GLOBIGNORE='*' :; name_array=($(cat $1))

i=0
array_len=${#name_array[@]}

while [ $i -lt $array_len ]
do
	bin/fastqtobam.py "${fastq_array[i]}" "${fastq_array[i+1]}" "${name_array[i]}" "${name_array[i+1]}"
	(( i=i+2 ))
done

IFS=$'\r\n' GLOBIGNORE='*' :; fs_bam_array=($(ls *.pe.srt.bam | cat))

i=0
while [ $i -lt $array_len  ]
do

	bin/human_contam.py -ap "${fs_bam_array[i]}" "${fs_bam_array[i+1]}"
	(( i=i+2 ))
done


today=$(date +"%m-%d-%Y")
mkdir "$today"
cd "$today"
mv ../*.bam* .
mv ../*.sam* .
mkdir basic_stats
mkdir control_stats
cp *.control*.filter.cs.bam control_stats
cp *.filter.cs.bam basic_stats

cd control_stats
IFS=$'\r\n' GLOBIGNORE='*' :; bam_array=($(ls *.bam | cat))
cp ../../*.genome .
number_control=${#bam_array[@]}

j=1
i=0
while [ $i -lt $number_control ]
do
	../../bin/control_stats.py "${bam_array[i]}" "${name_array[j]}.genome" > "${name_array[j-1]}.stats.txt"
	IFS=$'\r\n' GLOBIGNORE='*' :; num_chr=($(wc -l "${name_array[j]}.genome"))
	../../bin/control_graphs.R "${bam_array[i]:0:-4}" "$(( $num_chr-2))"
	(( i++ ))
	(( j=j+2 ))
done
mkdir graphs
mv *.pdf graphs

cd ../basic_stats
rm *.control*

IFS=$'\r\n' GLOBIGNORE='*' :; bam_array=($(ls | cat))
cp ../../*.genome .
cp ../../bed_reg_files/* .
IFS=$'\r\n' GLOBIGNORE='*' :; reg_array=($(ls *.bed* | cat))
array_len=${#bam_array[@]}

j=$(( $number_control*2+1))
i=0
while [ $i -lt $array_len ]
do
	../../bin/basic_statistic.py "${bam_array[i]}" "${reg_array[i]}" "${name_array[j]}.genome" > "${name_array[j-1]}.stats.txt"
	(( i++ ))
	(( j=j+2 ))
done


