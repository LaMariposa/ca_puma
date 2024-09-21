#!/bin/bash

#usage impute_chr.sh <chrID>

#script to impute a single chromosome
	#phase with shapeit
	#impute with glimpse

vcf=/scratch3/puma_megan/impute/infiles/Filt7_Post.vcf.gz
reflist=/scratch3/puma_megan/impute/infiles/ref155.list
targetlist=/scratch3/puma_megan/impute/infiles/lc329.list

threads=8

shapeit=/home/megan/bin/shapeit5/
glimpse=/home/megan/bin/GLIMPSE1/static_bins/

chr=$1


echo "analyzing chromosome $chr"


echo -e "\n\n\n ***** prepping reference panel *****"
bcftools view --samples-file ${reflist} -o refpanel_${chr}.vcf ${vcf} ${chr}
bgzip refpanel_${chr}.vcf
bcftools index refpanel_${chr}.vcf.gz


echo -e "\n\n\n ***** phasing reference panel with shapeit *****"
${shapeit}/phase_common_static --input refpanel_${chr}.vcf.gz --region ${chr} --output refpanel_${chr}.sphased.bcf --thread ${threads}
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' refpanel_${chr}.sphased.bcf | bgzip -c > refpanel_${chr}.sphased.tsv.gz
tabix -s1 -b2 -e2 refpanel_${chr}.sphased.tsv.gz 	


echo -e "\n\n\n ***** prep target panel *****"
bcftools view --samples-file ${targetlist} -o imppanel_${chr}.vcf ${vcf} ${chr}
bgzip imppanel_${chr}.vcf
bcftools index imppanel_${chr}.vcf.gz


echo -e "\n\n\n ***** imputing with glimpse using shapeit phased data *****"
${glimpse}/GLIMPSE_chunk_static --input refpanel_${chr}.sphased.bcf --region ${chr} --output schunks_${chr}.txt

while IFS="" read -r LINE || [ -n "$LINE" ];
do
printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
IRG=$(echo $LINE | cut -d" " -f3)
ORG=$(echo $LINE | cut -d" " -f4)
OUT=s_${chr}_${ID}.bcf
${glimpse}/GLIMPSE_phase_static --input imppanel_${chr}.vcf.gz --reference refpanel_${chr}.sphased.bcf --input-region ${IRG} --output-region ${ORG} --output ${OUT} --thread ${threads}	
bcftools index -f ${OUT}
done < schunks_${chr}.txt

ls s_${chr}_*.bcf > s_${chr}.list
${glimpse}/GLIMPSE_ligate_static --input s_${chr}.list --output ${chr}.sgimputed.bcf
bcftools index -f ${chr}.sgimputed.bcf
bcftools convert -O z -o ${chr}.sgimputed.vcf.gz ${chr}.sgimputed.bcf
tabix -p vcf ${chr}.sgimputed.vcf.gz

${glimpse}/GLIMPSE_sample_static --input ${chr}.sgimputed.bcf --solve --output ${chr}.sgimputedphased.bcf
bcftools index -f ${chr}.sgimputedphased.bcf


echo -e "\n\n\n ***** DONE!!! *****"
