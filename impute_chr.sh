#!/bin/bash

#usage impute_chr.sh <chrID>

#script to impute a single chromosome
	#phase with shapeit and beagle
	#impute (for both phasing) with glimpse

gmap=/scratch3/puma_megan/impute/inputfiles/puma_gentic.map
vcf=/scratch3/puma_megan/impute/inputfiles/downsampled_puma_SOI.vcf.gz
reflist=/scratch3/puma_megan/impute/inputfiles/testrefpanel125.list
targetlist=/scratch3/puma_megan/impute/inputfiles/lc329.list

threads=8

shapeit=/home/megan/bin/shapeit5/
beagle=/home/megan/bin/beagle.06Aug24.a91.jar
glimpse=/home/megan/bin/GLIMPSE1/static_bins/

chr=$1


echo "analyzing chromosome $chr"


echo -e "\n\n\n ***** prepping genetic maps *****"
grep ${chr} ${gmap} > puma_genetic_${chr}.map
awk '{print $2 "\t" $1 "\t" $3}' puma_genetic_${chr}.map > puma_poschcm.map
awk '{print $1 "\t.\t" $3 "\t" $2}' puma_genetic_${chr}.map | tail -n +2 > puma_plink.map


echo -e "\n\n\n ***** prepping reference panel *****"
bcftools view --samples-file ${reflist} -o refpanel_${chr}.vcf ${vcf} ${chr}
bgzip refpanel_${chr}.vcf
bcftools index refpanel_${chr}.vcf.gz


echo -e "\n\n\n ***** phasing reference panel with shapeit *****"
${shapeit}/phase_common_static --input refpanel_${chr}.vcf.gz --region ${chr} --map puma_poschcm.map --output refpanel_${chr}.sphased.bcf --thread ${threads}
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' refpanel_${chr}.sphased.bcf | bgzip -c > refpanel_${chr}.sphased.tsv.gz
tabix -s1 -b2 -e2 refpanel_${chr}.sphased.tsv.gz 	


echo -e "\n\n\n ***** phasing reference panel with beagle ***"
java -jar ${beagle} gt=refpanel_${chr}.vcf.gz map=puma_plink.map out=refpanel_${chr}.bphased nthreads=${threads}
	#could specify ne
bcftools index refpanel_${chr}.bphased.vcf.gz
bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' refpanel_${chr}.bphased.vcf.gz | bgzip -c > refpanel_${chr}.bphased.tsv.gz
tabix -s1 -b2 -e2 refpanel_${chr}.bphased.tsv.gz 

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
${glimpse}/GLIMPSE_phase_static --input imppanel_${chr}.vcf.gz --reference refpanel_${chr}.sphased.bcf --input-region ${IRG} \
	--map puma_poschcm.map \
	--output-region ${ORG} --output ${OUT} \
	--thread ${threads}
bcftools index -f ${OUT}
done < schunks_${chr}.txt

ls s_${chr}_*.bcf > s_${chr}.list
${glimpse}/GLIMPSE_ligate_static --input s_${chr}.list --output ${chr}.sgimputed.bcf
bcftools index -f ${chr}.sgimputed.bcf
bcftools convert -O z -o ${chr}.sgimputed.vcf.gz ${chr}.sgimputed.bcf
tabix -p vcf ${chr}.sgimputed.vcf.gz

${glimpse}/GLIMPSE_sample_static --input ${chr}.sgimputed.bcf --solve --output ${chr}.sgimputedphased.bcf
bcftools index -f ${chr}.sgimputedphased.bcf


echo -e "\n\n\n ***** imputing with glimpse using beagle phased data *****"
${glimpse}/GLIMPSE_chunk_static --input refpanel_${chr}.bphased.vcf.gz --region ${chr} --output bchunks_${chr}.txt

while IFS="" read -r LINE || [ -n "$LINE" ];
do
printf -v ID "%02d" $(echo $LINE | cut -d" " -f1)
IRG=$(echo $LINE | cut -d" " -f3)
ORG=$(echo $LINE | cut -d" " -f4)
OUT=b_${chr}_${ID}.bcf
${glimpse}/GLIMPSE_phase_static --input imppanel_${chr}.vcf.gz --reference refpanel_${chr}.bphased.vcf.gz --input-region ${IRG} \
	--map puma_poschcm.map \
	--output-region ${ORG} --output ${OUT} \
       	--thread ${threads}
bcftools index -f ${OUT}
done < bchunks_${chr}.txt

ls b_${chr}_*.bcf > b_${chr}.list
${glimpse}/GLIMPSE_ligate_static --input b_${chr}.list --output ${chr}.bgimputed.bcf
bcftools index -f ${chr}.bgimputed.bcf
bcftools convert -O z -o ${chr}.bgimputed.vcf.gz ${chr}.bgimputed.bcf
tabix -p vcf ${chr}.bgimputed.vcf.gz

${glimpse}/GLIMPSE_sample_static --input ${chr}.bgimputed.bcf --solve --output ${chr}.bgimputedphased.bcf
bcftools index -f ${chr}.bgimputedphased.bcf



echo -e "\n\n\n ***** DONE!!! *****"
