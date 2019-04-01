#-*- coding:UTF-8 -*-
import os
import sys
a=sys.argv[1]
os.system('date')
os.system('java -jar /home/sdb/biosoft/picardtools/CreateSequenceDictionary.jar R=reference_%s.fasta O=reference_%s.dict'%(a,a))
os.system('samtools merge filter_%s.bam bam_sorted_filter_A/*.bam'%(a))
os.system('samtools sort -o filter_sorted_%s.bam filter_%s.bam'%(a,a))
os.system('samtools index filter_sorted_%s.bam'%(a))
os.system('mkdir %s_HC_result'%(a))
os.system('java  -Xmx60G -jar /home/sdb/biosoft/gatk3.7/GenomeAnalysisTK.jar -T HaplotypeCaller --maxNumHaplotypesInPopulation 6 --dontUseSoftClippedBases -dt NONE --pcr_indel_model CONSERVATIVE -R reference_%s.fasta -I filter_sorted_%s.bam -o %s_HC_result/%s_gatk3.7_HC.vcf' %(a,a,a,a))
os.system('date')
