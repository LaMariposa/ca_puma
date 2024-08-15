#!/bin/bash

#usage impute_chr.sh <chrID>

#script to impute a single chromosome
	#phase with shapeit and beagle
	#impute (for both phasing) with glimpse

rmap=xxx
vcf=/scratch3/puma_megan/downsampled_puma_SOI.vcf.gz

reflist=testpanel.list
targetlist=lc329.list

threads=8


shapeit=/home/megan/bin/shapeit5/
glimpse=/home/megan/bin/GLIMPSE1/static_bins/

chr=$1


echo "analyzing chromosome $chr"


echo "prepping reference panel"
#bcftools view --samples-file ${reflist} -o refpanel_${chr}.vcf ${vcf} ${chr}
#bgzip refpanel_${chr}.vcf
#bcftools index refpanel_${chr}.vcf.gz


echo "phasing reference panel with shapeit"
#${shapeit}/phase_common_static --input refpanel_${chr}.vcf.gz --region ${chr} --output refpanel_${chr}.sphased.bcf --thread ${threads}
	#add --map ${rmap} when have rmap
#bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' refpanel_${chr}.sphased.bcf | bgzip -c > refpanel_${chr}.sphased.tsv.gz
#tabix -s1 -b2 -e2 refpanel_${chr}.sphased.tsv.gz 	


echo "prep target panel"
#bcftools view --samples-file lc329.list -o imppanel_${chr}.vcf ${vcf} ${chr}
#bgzip imppanel_${chr}.vcf
#bcftools index imppanel_${chr}.vcf.gz


echo "imputing shapeit phased with glimpse"
#${glimpse}/GLIMPSE_chunk_static --input refpanel_${chr}.sphased.bcf --region ${chr} --output schunks_${chr}.txt

while IFS="" read -r LINE || [ -n "$LINE" ];
do
printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
IRG=$(echo $LINE | cut -d" " -f3)
ORG=$(echo $LINE | cut -d" " -f4)
OUT=s_${chr}_${ID}.bcf
${glimpse}/GLIMPSE_phase_static --input imppanel_${chr}.vcf.gz --reference refpanel_${chr}.sphased.bcf --input-region ${IRG} \
	--output-region ${ORG} --output ${OUT} \
	--thread ${threads}
	#add --map chr20.b38.gmap.gz
bcftools index -f ${OUT}
done < schunks_${chr}.txt

