#!/bin/bash

# usage: rad_wgs_concord.sh <map> <radvcf> <wgsvcf>

#compare genotypes between rad and wgs


map=$1
radvcf=$2
wgsvcf=$3

while read -ra names
	do
		
		#get names
		radbam="${names[0]%\"}"
		wgsname="${names[1]%\"}"
		echo $wgsname
		radname="/pfs/tsfs1/project/wildgen/kgustaf7/combine/ddocent/WithErick/93bases/bam_pumcon1/${names[0]%\"}"
		echo $radname

		#compare genotypes
		java -jar /home/megan/bin/picard.jar GenotypeConcordance CALL_VCF=$radvcf CALL_SAMPLE=$radname O=$radbam TRUTH_VCF=$wgsvcf TRUTH_SAMPLE=$wgsname

       	done < $map
