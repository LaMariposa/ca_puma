#!/bin/bash

#usage impute_assess.sh <chr> <imputed.vcf>

#script to assess imputation using downsampled samples and rad data

chr=$1
impvcf=$2

refvcf=/scratch3/puma_megan/impute/assess/refpanel_${chr}.vcf.gz
goldvcf=/scratch3/puma_megan/impute/assess/gold_${chr}.vcf.gz
gold_ds_names=/scratch3/puma_megan/impute/inputfiles/matchdsnames.txt
dssamples=/scratch3/puma_megan/impute/inputfiles/downsample30.list

threads=8

glimpse=/home/megan/bin/GLIMPSE1/static_bins/



echo "analyzing chromosome $chr and imputed vcf $impvcf"


echo -e "\n\n\n ***** prep input files *****"
bcftools reheader --samples ${gold_ds_names} -o imp_renamed.vcf.gz ${impvcf}
bcftools index imp_renamed.vcf.gz
echo ${chr} ${refvcf} ${goldvcf} imp_renamed.vcf.gz > infiles.txt


echo -e "\n\n\n ***** assess concordance *****"
${glimpse}/GLIMPSE_concordance_static --input infiles.txt --samples ${dssamples} --gt-validation --gt-target --minPROB 0.9999 --minDP 0 --output concord_out --thread ${threads}  --bins 0.00000 0.00100 0.00200 0.00500 0.01000 0.05000 0.10000 0.20000 0.50000

/scratch3/puma_megan/impute/assess/concordance_plot.py
prefix=${PWD##*/}
mv accplot.png ${prefix}_concord.png

