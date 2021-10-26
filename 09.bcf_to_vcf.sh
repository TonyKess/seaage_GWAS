module load bcftools

cd scratch/Salmon_AOM/seaage/CIGENE_align/

 for i in {01..29} ; do bcftools view Salmo_ssa$i\_CIGENE_80geno.bcf > Salmo_ssa$i\_CIGENE_80geno.vcf; done


