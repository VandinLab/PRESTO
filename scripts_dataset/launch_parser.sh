#!/bin/bash 

FOLDER=$1
DEST="first_pass/"
a=1
for capture in `ls ${PWD}/${FOLDER} | sort -V`
do
	python parse_pcap.py "${FOLDER}${capture}"  > ${DEST}file_${a}
	((a++))
done
