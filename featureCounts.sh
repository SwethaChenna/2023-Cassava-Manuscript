#### De novo annotation

gff=./Called_transcripts_both_WTonly_Cold3hOnly_sorted_mod.gff


#### ssRNA-Seq
featureCounts -p -s 2 -M -T 30 --countReadPairs -g Parent -t exon -a $gff -o WT_Cold3h_ssRNA-Seq_featureCounts.txt ../ssRNA-Seq/v6.54/Bam/sR1_WT_Rep1_Aligned.out_mapq_fixmate_sorted_dedup.bam ../ssRNA-Seq/v6.54/Bam/sR2_WT_Rep2_Aligned.out_mapq_fixmate_sorted_dedup.bam ../ssRNA-Seq/v6.54/Bam/sR3_WT_Rep3_Aligned.out_mapq_fixmate_sorted_dedup.bam ../ssRNA-Seq/v6.54/Bam/sR4_WT_Rep4_Aligned.out_mapq_fixmate_sorted_dedup.bam ../ssRNA-Seq/v6.54/Bam/sR5_Cold3h_Rep1_Aligned.out_mapq_fixmate_sorted_dedup.bam ../ssRNA-Seq/v6.54/Bam/sR6_Cold3h_Rep2_Aligned.out_mapq_fixmate_sorted_dedup.bam ../ssRNA-Seq/v6.54/Bam/sR7_Cold3h_Rep3_Aligned.out_mapq_fixmate_sorted_dedup.bam ../ssRNA-Seq/v6.54/Bam/sR8_Cold3h_Rep4_Aligned.out_mapq_fixmate_sorted_dedup.bam


#### Quant-Seq
featureCounts -s 2 -M -T 30 -g Parent -t exon -a $gff -o WT_Cold3h_Quant-Seq_featureCounts.txt ../Quant-Seq/v6.54/Bam/Cold3h_rep1_Aligned.out_mapq.bam ../Quant-Seq/v6.54/Bam/Cold3h_rep2_Aligned.out_mapq.bam ../Quant-Seq/v6.54/Bam/Cold3h_rep3_Aligned.out_mapq.bam ../Quant-Seq/v6.54/Bam/Cold3h_rep4_Aligned.out_mapq.bam ../Quant-Seq/v6.54/Bam/TM12_rep1_Aligned.out_mapq.bam ../Quant-Seq/v6.54/Bam/TM12_rep2_Aligned.out_mapq.bam ../Quant-Seq/v6.54/Bam/TM12_rep3_Aligned.out_mapq.bam ../Quant-Seq/v6.54/Bam/TM12_rep4_Aligned.out_mapq.bam


#### STRIPE-Seq
featureCounts -s 1 -M -T 30 -g Parent -t exon -a $gff -o WT_Cold3h_STRIPE-Seq_featureCounts.txt ../STRIPE-Seq/v6.54/Bam/TM12_rep1_UMI_trimmed_Aligned.out_mapq_dedup_sorted.bam ../STRIPE-Seq/v6.54/Bam/TM12_rep2_UMI_trimmed_Aligned.out_mapq_dedup_sorted.bam ../STRIPE-Seq/v6.54/Bam/TM12_rep3_UMI_trimmed_Aligned.out_mapq_dedup_sorted.bam ../STRIPE-Seq/v6.54/Bam/TM12_rep4_UMI_trimmed_Aligned.out_mapq_dedup_sorted.bam ../STRIPE-Seq/v6.54/Bam/Cold3h_rep1_UMI_trimmed_Aligned.out_mapq_dedup_sorted.bam ../STRIPE-Seq/v6.54/Bam/Cold3h_rep2_UMI_trimmed_Aligned.out_mapq_dedup_sorted.bam ../STRIPE-Seq/v6.54/Bam/Cold3h_rep3_UMI_trimmed_Aligned.out_mapq_dedup_sorted.bam ../STRIPE-Seq/v6.54/Bam/Cold3h_rep4_UMI_trimmed_Aligned.out_mapq_dedup_sorted.bam



#### EnsemBlPlants

gff=../../Genomes/Cassava/v6.54/Manihot_esculenta_v6.54_woScaffold_filtered.gff3


#### ssRNA-Seq
featureCounts -p -s 2 -M -T 30 --countReadPairs -g Parent -t exon -a $gff -o WT_Cold3h_ssRNA-Seq_featureCounts_Ensembl.txt ../ssRNA-Seq/v6.54/Bam/sR1_WT_Rep1_Aligned.out_mapq_fixmate_sorted_dedup.bam ../ssRNA-Seq/v6.54/Bam/sR2_WT_Rep2_Aligned.out_mapq_fixmate_sorted_dedup.bam ../ssRNA-Seq/v6.54/Bam/sR3_WT_Rep3_Aligned.out_mapq_fixmate_sorted_dedup.bam ../ssRNA-Seq/v6.54/Bam/sR4_WT_Rep4_Aligned.out_mapq_fixmate_sorted_dedup.bam ../ssRNA-Seq/v6.54/Bam/sR5_Cold3h_Rep1_Aligned.out_mapq_fixmate_sorted_dedup.bam ../ssRNA-Seq/v6.54/Bam/sR6_Cold3h_Rep2_Aligned.out_mapq_fixmate_sorted_dedup.bam ../ssRNA-Seq/v6.54/Bam/sR7_Cold3h_Rep3_Aligned.out_mapq_fixmate_sorted_dedup.bam ../ssRNA-Seq/v6.54/Bam/sR8_Cold3h_Rep4_Aligned.out_mapq_fixmate_sorted_dedup.bam


#### Quant-Seq
featureCounts -s 2 -M -T 30 -g Parent -t exon -a $gff -o WT_Cold3h_Quant-Seq_featureCounts_Ensembl.txt ../Quant-Seq/v6.54/Bam/Cold3h_rep1_Aligned.out_mapq.bam ../Quant-Seq/v6.54/Bam/Cold3h_rep2_Aligned.out_mapq.bam ../Quant-Seq/v6.54/Bam/Cold3h_rep3_Aligned.out_mapq.bam ../Quant-Seq/v6.54/Bam/Cold3h_rep4_Aligned.out_mapq.bam ../Quant-Seq/v6.54/Bam/TM12_rep1_Aligned.out_mapq.bam ../Quant-Seq/v6.54/Bam/TM12_rep2_Aligned.out_mapq.bam ../Quant-Seq/v6.54/Bam/TM12_rep3_Aligned.out_mapq.bam ../Quant-Seq/v6.54/Bam/TM12_rep4_Aligned.out_mapq.bam


#### STRIPE-Seq
featureCounts -s 1 -M -T 30 -g Parent -t exon -a $gff -o WT_Cold3h_STRIPE-Seq_featureCounts_Ensembl.txt ../STRIPE-Seq/v6.54/Bam/TM12_rep1_UMI_trimmed_Aligned.out_mapq_dedup_sorted.bam ../STRIPE-Seq/v6.54/Bam/TM12_rep2_UMI_trimmed_Aligned.out_mapq_dedup_sorted.bam ../STRIPE-Seq/v6.54/Bam/TM12_rep3_UMI_trimmed_Aligned.out_mapq_dedup_sorted.bam ../STRIPE-Seq/v6.54/Bam/TM12_rep4_UMI_trimmed_Aligned.out_mapq_dedup_sorted.bam ../STRIPE-Seq/v6.54/Bam/Cold3h_rep1_UMI_trimmed_Aligned.out_mapq_dedup_sorted.bam ../STRIPE-Seq/v6.54/Bam/Cold3h_rep2_UMI_trimmed_Aligned.out_mapq_dedup_sorted.bam ../STRIPE-Seq/v6.54/Bam/Cold3h_rep3_UMI_trimmed_Aligned.out_mapq_dedup_sorted.bam ../STRIPE-Seq/v6.54/Bam/Cold3h_rep4_UMI_trimmed_Aligned.out_mapq_dedup_sorted.bam
