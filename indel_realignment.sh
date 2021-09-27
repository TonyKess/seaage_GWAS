## Load software modules
module load nixpkgs/16.09
module load gatk/3.7

#First index dedup reads - do with Salmo_AOM_dedupindex.sh
cd /home/tkess/scratch/Salmon_AOM/seaage/CIGENE_align/

#get to realigning
cat ../set12.txt | \
  parallel --jobs 10 \
 'java -Xmx10G -jar $EBROOTGATK/GenomeAnalysisTK.jar \
 -T IndelRealigner \
 -R /home/tkess/scratch/Salmon_CIGENE_Ref/reference/CIGENE-ICSASG_v2.fa \
 -I {}.deDup.bam  \
 -targetIntervals {}.intervals \
 -o {}.realigned.bam '
