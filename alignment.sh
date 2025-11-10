#!/bin/bash -l        

# Alignment for P24X0

STAR --runThreadN 12 --runMode genomeGenerate --genomeDir ref/P24XO_Male --genomeFastaFiles /ref/P24XO_Male/Vv_P24XO_Male_ref.fasta --sjdbGTFfile ref/P24XO_Male/P24XO_Male.gff --sjdbOverhang 1 --sjdbGTFfeatureExon CDS

# Male Testis
for i in S{53..60}; 
	do 
	STAR --runMode alignReads --outSAMtype BAM SortedByCoordinate --runThreadN 12 --readFilesCommand zcat --genomeDir ref/P24XO_Male/ --sjdbGTFfile ref/P24XO_Male/P24XO_Male.gff --sjdbGTFfeatureExon CDS --quantMode TranscriptomeSAM GeneCounts --bamRemoveDuplicatesType UniqueIdenticalNotMulti --outFileNamePrefix Male_testis_"$i"_onP24XOM --readFilesIn P24XO_samples/male_testis/"$i"_R1.fastq.gz P24XO_samples/male_testis/"$i"_R2.fastq.gz;
	bamtools index -in Male_testis_"$i"_onP24XOMAligned.sortedByCoord.out.bam;
done
featureCounts -p -B -T 5 -F GFF -g ID -t gene -a ref/P24XO_Male/P24XO_Male.gff -o Male_testis_all_onP24XOM_featureCounts_new Male_testis_S53_onP24XOMAligned.sortedByCoord.out.bam Male_testis_S54_onP24XOMAligned.sortedByCoord.out.bam Male_testis_S55_onP24XOMAligned.sortedByCoord.out.bam Male_testis_S56_onP24XOMAligned.sortedByCoord.out.bam Male_testis_S57_onP24XOMAligned.sortedByCoord.out.bam Male_testis_S58_onP24XOMAligned.sortedByCoord.out.bam Male_testis_S59_onP24XOMAligned.sortedByCoord.out.bam Male_testis_S60_onP24XOMAligned.sortedByCoord.out.bam

# Female ovaries
for i in S{21..30};
	do 
	STAR --runMode alignReads --outSAMtype BAM SortedByCoordinate --runThreadN 12 --readFilesCommand zcat --genomeDir /crex/proj/snic2021-6-151/private/Suvi/ref/P24XO_Male/ --sjdbGTFfile /crex/proj/snic2021-6-151/private/Suvi/ref/P24XO_Male/P24XO_Male.gff --sjdbGTFfeatureExon CDS --quantMode TranscriptomeSAM GeneCounts --bamRemoveDuplicatesType UniqueIdenticalNotMulti --outFileNamePrefix Female_ovaries_"$i"_onP24XOM --readFilesIn P24XO_samples/female_ovaries/"$i"_R1.fastq.gz P24XO_samples/female_ovaries/"$i"_R2.fastq.gz;
	bamtools index -in Female_ovaries_"$i"_onP24XOMAligned.sortedByCoord.out.bam;
done
featureCounts -p -B -T 5 -F GFF -g ID -t gene -a ref/P24XO_Male/P24XO_Male.gff -o Female_ovaries_all_onP24XOM_featureCounts_new Female_ovaries_S21_onP24XOMAligned.sortedByCoord.out.bam Female_ovaries_S22_onP24XOMAligned.sortedByCoord.out.bam Female_ovaries_S23_onP24XOMAligned.sortedByCoord.out.bam Female_ovaries_S24_onP24XOMAligned.sortedByCoord.out.bam Female_ovaries_S25_onP24XOMAligned.sortedByCoord.out.bam Female_ovaries_S26_onP24XOMAligned.sortedByCoord.out.bam Female_ovaries_S27_onP24XOMAligned.sortedByCoord.out.bam Female_ovaries_S28_onP24XOMAligned.sortedByCoord.out.bam Female_ovaries_S30_onP24XOMAligned.sortedByCoord.out.bam Female_ovaries_S30_onP24XOMAligned.sortedByCoord.out.bam

# Female head 
for i in S{1..10};
	do 
	STAR --runMode alignReads --outSAMtype BAM SortedByCoordinate --runThreadN 12 --readFilesCommand zcat --genomeDir ref/P24XO_Male/ --sjdbGTFfile /ref/P24XO_Male/P24XO_Male.gff --sjdbGTFfeatureExon CDS --quantMode TranscriptomeSAM GeneCounts --bamRemoveDuplicatesType UniqueIdenticalNotMulti --outFileNamePrefix Female_head_"$i"_onP24XOM --readFilesIn P24XO_samples/female_head/"$i"_R1.fastq.gz P24XO_samples/female_head/"$i"_R2.fastq.gz;
	bamtools index -in Female_head_"$i"_onP24XOMAligned.sortedByCoord.out.bam;
done
featureCounts -p -B -T 5 -F GFF -g ID -t gene -a ref/P24XO_Male/P24XO_Male.gff -o Female_head_all_onP24XOM_featureCounts_new Female_head_S1_onP24XOMAligned.sortedByCoord.out.bam Female_head_S2_onP24XOMAligned.sortedByCoord.out.bam Female_head_S3_onP24XOMAligned.sortedByCoord.out.bam Female_head_S4_onP24XOMAligned.sortedByCoord.out.bam Female_head_S5_onP24XOMAligned.sortedByCoord.out.bam Female_head_S6_onP24XOMAligned.sortedByCoord.out.bam Female_head_S7_onP24XOMAligned.sortedByCoord.out.bam Female_head_S8_onP24XOMAligned.sortedByCoord.out.bam Female_head_S9_onP24XOMAligned.sortedByCoord.out.bam Female_head_S10_onP24XOMAligned.sortedByCoord.out.bam

# Male head
for i in S{31..38};
	do 
	STAR --runMode alignReads --outSAMtype BAM SortedByCoordinate --runThreadN 12 --readFilesCommand zcat --genomeDir ref/P24XO_Male/ --sjdbGTFfile ref/P24XO_Male/P24XO_Male.gff --sjdbGTFfeatureExon CDS --quantMode TranscriptomeSAM GeneCounts --bamRemoveDuplicatesType UniqueIdenticalNotMulti --outFileNamePrefix Male_head_"$i"_onP24XOM --readFilesIn P24XO_samples/male_head/"$i"_R1.fastq.gz P24XO_samples/male_head/"$i"_R2.fastq.gz;
	bamtools index -in Male_head_"$i"_onP24XOMAligned.sortedByCoord.out.bam;
done
featureCounts -p -B -T 5 -F GFF -g ID -t gene -a ref/P24XO_Male/P24XO_Male.gff -o Male_head_all_onP24XOM_featureCounts_new Male_head_S31_onP24XOMAligned.sortedByCoord.out.bam Male_head_S32_onP24XOMAligned.sortedByCoord.out.bam Male_head_S33_onP24XOMAligned.sortedByCoord.out.bam Male_head_S34_onP24XOMAligned.sortedByCoord.out.bam Male_head_S35_onP24XOMAligned.sortedByCoord.out.bam Male_head_S36_onP24XOMAligned.sortedByCoord.out.bam Male_head_S37_onP24XOMAligned.sortedByCoord.out.bam Male_head_S38_onP24XOMAligned.sortedByCoord.out.bam

# Female legs
for i in S{11..20};
	do 
	STAR --runMode alignReads --outSAMtype BAM SortedByCoordinate --runThreadN 12 --readFilesCommand zcat --genomeDir ref/P24XO_Male/ --sjdbGTFfile ref/P24XO_Male/P24XO_Male.gff --sjdbGTFfeatureExon CDS --quantMode TranscriptomeSAM GeneCounts --bamRemoveDuplicatesType UniqueIdenticalNotMulti --outFileNamePrefix Female_legs_"$i"_onP24XOM --readFilesIn P24XO_samples/female_legs/"$i"_R1.fastq.gz P24XO_samples/female_legs/"$i"_R2.fastq.gz;
	bamtools index -in Female_legs_"$i"_onP24XOMAligned.sortedByCoord.out.bam;
done
featureCounts -p -B -T 5 -F GFF -g ID -t gene -a ref/P24XO_Male/P24XO_Male.gff -o Female_legs_all_onP24XOM_featureCounts_new Female_legs_S11_onP24XOMAligned.sortedByCoord.out.bam Female_legs_S12_onP24XOMAligned.sortedByCoord.out.bam Female_legs_S13_onP24XOMAligned.sortedByCoord.out.bam Female_legs_S14_onP24XOMAligned.sortedByCoord.out.bam Female_legs_S15_onP24XOMAligned.sortedByCoord.out.bam Female_legs_S16_onP24XOMAligned.sortedByCoord.out.bam Female_legs_S17_onP24XOMAligned.sortedByCoord.out.bam Female_legs_S18_onP24XOMAligned.sortedByCoord.out.bam Female_legs_S19_onP24XOMAligned.sortedByCoord.out.bam Female_legs_S20_onP24XOMAligned.sortedByCoord.out.bam

# Male legs
for i in S{39..45};
	do 
	STAR --runMode alignReads --outSAMtype BAM SortedByCoordinate --runThreadN 12 --readFilesCommand zcat --genomeDir ref/P24XO_Male/ --sjdbGTFfile ref/P24XO_Male/P24XO_Male.gff --sjdbGTFfeatureExon CDS --quantMode TranscriptomeSAM GeneCounts --bamRemoveDuplicatesType UniqueIdenticalNotMulti --outFileNamePrefix Male_legs_"$i"_onP24XOM --readFilesIn P24XO_samples/male_legs/"$i"_R1.fastq.gz P24XO_samples/male_legs/"$i"_R2.fastq.gz;
	bamtools index -in Male_legs_"$i"_onP24XOMAligned.sortedByCoord.out.bam;
done
featureCounts -p -B -T 5 -F GFF -g ID -t gene -a /ref/P24XO_Male/P24XO_Male.gff -o Male_legs_all_onP24XOM_featureCounts_new Male_legs_S39_onP24XOMAligned.sortedByCoord.out.bam Male_legs_S40_onP24XOMAligned.sortedByCoord.out.bam Male_legs_S41_onP24XOMAligned.sortedByCoord.out.bam Male_legs_S42_onP24XOMAligned.sortedByCoord.out.bam Male_legs_S43_onP24XOMAligned.sortedByCoord.out.bam Male_legs_S44_onP24XOMAligned.sortedByCoord.out.bam Male_legs_S45_onP24XOMAligned.sortedByCoord.out.bam

# Alignment for P24XY

STAR --runThreadN 12 --runMode genomeGenerate --genomeDir ref/P24XY_Female --genomeFastaFiles ref/P24XY_Female/Vv_P24XY_Female_ref.fasta --sjdbGTFfile ref/P24XY_Female/P24XY_Female.gff --sjdbOverhang 1 --sjdbGTFfeatureExon CDS

# Male Testis 
for i in S{98..104}; 
	do 
	STAR --runMode alignReads --outSAMtype BAM SortedByCoordinate --runThreadN 12 --readFilesCommand zcat --genomeDir ref/P24XY_Female/ --sjdbGTFfile ref/P24XY_Female/P24XY_Female.gff --sjdbGTFfeatureExon CDS --quantMode TranscriptomeSAM GeneCounts --bamRemoveDuplicatesType UniqueIdenticalNotMulti --outFileNamePrefix Male_testis_"$i"_onP24XYF --readFilesIn P24XY_samples/male_testis/"$i"_R1.fastq.gz P24XY_samples/male_testis/"$i"_R2.fastq.gz;
	bamtools index -in Male_testis_"$i"_onP24XYFAligned.sortedByCoord.out.bam
done
featureCounts -p -B -T 5 -F GFF -g ID -s 2 -t gene -a ref/P24XY_Female/P24XY_Female.gff -o Male_testis_all_onP24XYF_featureCounts Male_testis_S98_onP24XYFAligned.sortedByCoord.out.bam Male_testis_S99_onP24XYFAligned.sortedByCoord.out.bam Male_testis_S100_onP24XYFAligned.sortedByCoord.out.bam Male_testis_S101_onP24XYFAligned.sortedByCoord.out.bam Male_testis_S102_onP24XYFAligned.sortedByCoord.out.bam Male_testis_S103_onP24XYFAligned.sortedByCoord.out.bam Male_testis_S104_onP24XYFAligned.sortedByCoord.out.bam

# Female ovaries
for i in S{51..52} S{77..85}
	do
  STAR --runMode alignReads --outSAMtype BAM SortedByCoordinate --runThreadN 12 --readFilesCommand zcat --genomeDir ref/P24XY_Female/ --sjdbGTFfile ref/P24XY_Female/P24XY_Female.gff --sjdbGTFfeatureExon CDS --quantMode TranscriptomeSAM GeneCounts --bamRemoveDuplicatesType UniqueIdenticalNotMulti --outFileNamePrefix Female_ovaries_"$i"_onP24XYF --readFilesIn P24XY_samples/female_ovaries/"$i"_R1.fastq.gz P24XY_samples/female_ovaries/"$i"_R2.fastq.gz;
	bamtools index -in Female_ovaries_"$i"_onP24XYFAligned.sortedByCoord.out.bam
done
featureCounts -p -B -T 5 -F GFF -g ID -s 2 -t gene -a ref/P24XY_Female/P24XY_Female.gff -o Female_ovaries_all_onP24XYF_featureCounts Female_ovaries_S51_onP24XYFAligned.sortedByCoord.out.bam Female_ovaries_S52_onP24XYFAligned.sortedByCoord.out.bam Female_ovaries_S77_onP24XYFAligned.sortedByCoord.out.bam Female_ovaries_S78_onP24XYFAligned.sortedByCoord.out.bam Female_ovaries_S79_onP24XYFAligned.sortedByCoord.out.bam Female_ovaries_S80_onP24XYFAligned.sortedByCoord.out.bam Female_ovaries_S81_onP24XYFAligned.sortedByCoord.out.bam Female_ovaries_S82_onP24XYFAligned.sortedByCoord.out.bam Female_ovaries_S83_onP24XYFAligned.sortedByCoord.out.bam Female_ovaries_S84_onP24XYFAligned.sortedByCoord.out.bam Female_ovaries_S85_onP24XYFAligned.sortedByCoord.out.bam

# Female head
for i in S{46..50} S{61..66};
	do 
  STAR --runMode alignReads --outSAMtype BAM SortedByCoordinate --runThreadN 12 --readFilesCommand zcat --genomeDir ref/P24XY_Female/ --sjdbGTFfile ref/P24XY_Female/P24XY_Female.gff --sjdbGTFfeatureExon CDS --quantMode TranscriptomeSAM GeneCounts --bamRemoveDuplicatesType UniqueIdenticalNotMulti --outFileNamePrefix Female_head_"$i"_onP24XYF --readFilesIn P24XY_samples/female_head/"$i"_R1.fastq.gz P24XY_samples/female_head/"$i"_R2.fastq.gz;
	bamtools index -in Female_head_"$i"_onP24XYFAligned.sortedByCoord.out.bam
  done
  featureCounts -p -B -T 5 -F GFF -g ID -s 2 -t gene -a ref/P24XY_Female/P24XY_Female.gff -o Female_head_all_onP24XYF_featureCounts Female_head_S46_onP24XYFAligned.sortedByCoord.out.bam Female_head_S47_onP24XYFAligned.sortedByCoord.out.bam Female_head_S48_onP24XYFAligned.sortedByCoord.out.bam Female_head_S49_onP24XYFAligned.sortedByCoord.out.bam Female_head_S50_onP24XYFAligned.sortedByCoord.out.bam Female_head_S61_onP24XYFAligned.sortedByCoord.out.bam Female_head_S62_onP24XYFAligned.sortedByCoord.out.bam Female_head_S63_onP24XYFAligned.sortedByCoord.out.bam Female_head_S64_onP24XYFAligned.sortedByCoord.out.bam Female_head_S65_onP24XYFAligned.sortedByCoord.out.bam Female_head_S66_onP24XYFAligned.sortedByCoord.out.bam

# Male head
for i in S{86..92};
	do
	STAR --runMode alignReads --outSAMtype BAM SortedByCoordinate --runThreadN 12 --readFilesCommand zcat --genomeDir ref/P24XY_Female/ --sjdbGTFfile /ref/P24XY_Female/P24XY_Female.gff --sjdbGTFfeatureExon CDS --quantMode TranscriptomeSAM GeneCounts --bamRemoveDuplicatesType UniqueIdenticalNotMulti --outFileNamePrefix Male_head_"$i"_onP24XYF --readFilesIn P24XY_samples/male_head/"$i"_R1.fastq.gz P24XY_samples/male_head/"$i"_R2.fastq.gz;
	bamtools index -in Male_head_"$i"_onP24XYFAligned.sortedByCoord.out.bam
done
featureCounts -p -B -T 5 -F GFF -g ID -s 2 -t gene -a ref/P24XY_Female/P24XY_Female.gff -o Male_head_all_onP24XYF_featureCounts Male_head_S86_onP24XYFAligned.sortedByCoord.out.bam Male_head_S87_onP24XYFAligned.sortedByCoord.out.bam Male_head_S88_onP24XYFAligned.sortedByCoord.out.bam Male_head_S89_onP24XYFAligned.sortedByCoord.out.bam Male_head_S90_onP24XYFAligned.sortedByCoord.out.bam Male_head_S91_onP24XYFAligned.sortedByCoord.out.bam Male_head_S92_onP24XYFAligned.sortedByCoord.out.bam

# Female legs
for i in S{67..76};
	do 
	STAR --runMode alignReads --outSAMtype BAM SortedByCoordinate --runThreadN 12 --readFilesCommand zcat --genomeDir ref/P24XY_Female/ --sjdbGTFfile ref/P24XY_Female/P24XY_Female.gff --sjdbGTFfeatureExon CDS --quantMode TranscriptomeSAM GeneCounts --bamRemoveDuplicatesType UniqueIdenticalNotMulti --outFileNamePrefix Female_legs_"$i"_onP24XYF --readFilesIn P24XY_samples/female_legs/"$i"_R1.fastq.gz P24XY_samples/female_legs/"$i"_R2.fastq.gz;
	bamtools index -in Female_legs_"$i"_onP24XYFAligned.sortedByCoord.out.bam;
done
featureCounts -p -B -T 5 -F GFF -g ID -s 2 -t gene -a ref/P24XY_Female/P24XY_Female.gff -o Female_legs_all_onP24XYF_featureCounts Female_legs_S67_onP24XYFAligned.sortedByCoord.out.bam Female_legs_S68_onP24XYFAligned.sortedByCoord.out.bam Female_legs_S69_onP24XYFAligned.sortedByCoord.out.bam Female_legs_S70_onP24XYFAligned.sortedByCoord.out.bam Female_legs_S71_onP24XYFAligned.sortedByCoord.out.bam Female_legs_S72_onP24XYFAligned.sortedByCoord.out.bam Female_legs_S73_onP24XYFAligned.sortedByCoord.out.bam Female_legs_S74_onP24XYFAligned.sortedByCoord.out.bam Female_legs_S75_onP24XYFAligned.sortedByCoord.out.bam Female_legs_S76_onP24XYFAligned.sortedByCoord.out.bam

# Male legs 
for i in S{93..97};
	do 
	STAR --runMode alignReads --outSAMtype BAM SortedByCoordinate --runThreadN 12 --readFilesCommand zcat --genomeDir ref/P24XY_Female/ --sjdbGTFfile ref/P24XY_Female/P24XY_Female.gff --sjdbGTFfeatureExon CDS --quantMode TranscriptomeSAM GeneCounts --bamRemoveDuplicatesType UniqueIdenticalNotMulti --outFileNamePrefix Male_legs_"$i"_onP24XYF --readFilesIn P24XY_samples/male_legs/"$i"_R1.fastq.gz P24XY_samples/male_legs/"$i"_R2.fastq.gz;
	bamtools index -in Male_legs_"$i"_onP24XYFAligned.sortedByCoord.out.bam;
done
featureCounts -p -B -T 5 -F GFF -g ID -s 2 -t gene -a ref/P24XY_Female/P24XY_Female.gff -o Male_legs_all_onP24XYF_featureCounts Male_legs_S93_onP24XYFAligned.sortedByCoord.out.bam Male_legs_S94_onP24XYFAligned.sortedByCoord.out.bam Male_legs_S95_onP24XYFAligned.sortedByCoord.out.bam Male_legs_S96_onP24XYFAligned.sortedByCoord.out.bam Male_legs_S97_onP24XYFAligned.sortedByCoord.out.bam
