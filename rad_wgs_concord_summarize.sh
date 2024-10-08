#!/bin/bash

#compare genotypes between rad and wgs

# usage  rad_wgs_concord_summarize.sh <id map>

map=$1
#echo $map
#map=rad_wgs_names_comps1.tsv

echo -e "ID\tnmatches" > matches.txt
echo -e "ID\tnmismatches\tconcordance" > mismatches.txt

while read -ra names
	do
		
		#get names
		wgsname="${names[1]%\"}"
		#echo $wgsname
		radname="${names[0]%\"}"
		#echo $radname

		#get matches
		grep "SNP" ${radname}.genotype_concordance_detail_metrics | grep -v "MISSING" | grep -v "NO_CALL" | grep -v "VC_FILTERED" | \
			awk '{if ($4==$5) print $0}' | awk -v OFS='\t' '{pos+=$6} END{print $2,pos}' >> matches.txt
       
                #get mismatches
                grep "SNP" ${radname}.genotype_concordance_detail_metrics | grep -v "MISSING" | grep -v "NO_CALL" | grep -v "VC_FILTERED" | \
                        awk '{if ($4!=$5) print $0}' | awk -v OFS='\t' '{neg+=$6} END{print $2, neg}' >> mismatches.txt

	done < $map

	#combine and calc final
	paste matches.txt mismatches.txt | awk -v OFS='\t' '{tot=$2+$4; conc=NA; if(tot>0){conc=$2/tot}; print $0,conc}'
