# Step 3: RNA-seq Variant Calling Workflow (GATK)

Follows the [GATK Best Practices workflow for RNA-seq short variant discovery](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-),
starting from STAR-aligned BAM files. Tools used: `samtools`, `GATK 4.3.0.0`,
`picard`, `snpEff`, and `BEDOPS`.

> Two sample groups are processed with the same steps below — a single-sample run
> (e.g. sample `94`) and a merged-replicate run (e.g. samples `94`/`95`/`96` merged, and
> `97`/`98`/`99` merged). Replace file names/paths with your own sample IDs.

## 1. Merge BAM files from replicates

Merge BAM files from three replicates that were generated after STAR alignment.

```bash
samtools merge MI949596sorted_out.bam ~/Desktop/MIA1bamfiles/*.bam
```

## 2. Mark duplicates

```bash
gatk MarkDuplicates \
     -I Aligned94.sortedByCoord.out.bam \
     -O Aligned94.dedup.bam \
     -M Aligned94.dedup.txt
```

## 3. Add read group (RG) tags

```bash
java -jar picard.jar AddOrReplaceReadGroups \
     I=Aligned94.dedup.bam \
     O=Aligned94.dedup.RGI.bam \
     RGID=1 RGLB=lib2 RGPL=illumina RGPU=unit1 RGSM=3
```

## 4. Create the reference sequence dictionary

```bash
gatk CreateSequenceDictionary \
     -R Microbotryum_intermedium_gca_900096595.Microbotrium_Intermedium_Assembly.dna.toplevel.fa
```

## 5. Index the reference sequence

```bash
samtools faidx Microbotryum_intermedium_gca_900096595.Microbotrium_Intermedium_Assembly.dna.toplevel.fa
```

## 6. Handle splicing events (SplitNCigarReads)

```bash
gatk SplitNCigarReads \
     -R Microbotryum_intermedium_gca_900096595.Microbotrium_Intermedium_Assembly.dna.toplevel.fa \
     -I Aligned94.dedup.RGI.bam \
     -O Aligned94.dedup_reads_RGI-splitreads.bam
```

## 7. Variant calling (HaplotypeCaller)

```bash
gatk --java-options "-Xmx4g" HaplotypeCaller \
     -R Microbotryum_intermedium_gca_900096595.Microbotrium_Intermedium_Assembly.dna.toplevel.fa \
     -I Aligned94.dedup_reads_RGI-splitreads.bam \
     -O raw_varients.vcf \
     -bamout bamout.bam
```

---

### Merged-replicate run (steps 8–14)

## 8. Merge BAM files (replicate set)

```bash
samtools merge MI949596sorted_out.bam ~/Desktop/MIA1bamfiles/*.bam
```

## 9. Mark duplicates

```bash
gatk MarkDuplicates \
     -I MI949596sorted_out.bam \
     -O MI949596sorted_out.dedup.bam \
     -M Metrix949596.dedup.txt
```

## 10. Add read group (RG) tags

```bash
java -jar picard.jar AddOrReplaceReadGroups \
     I=MI949596sorted_out.dedup.bam \
     O=MI949596sorted_dedup_RGI.bam \
     RGID=1 RGLB=lib2 RGPL=illumina RGPU=unit1 RGSM=3
```

## 11. Handle splicing events (SplitNCigarReads)

```bash
gatk SplitNCigarReads \
     -R Microbotryum_intermedium_gca_900096595.Microbotrium_Intermedium_Assembly.dna.toplevel.fa \
     -I MI949596sorted_dedup_RGI.bam \
     -O MI949596sorted_dedup_RGI-splitreads.bam
```

## 12. Variant calling (HaplotypeCaller)

```bash
gatk --java-options "-Xmx4g" HaplotypeCaller \
     -R Microbotryum_intermedium_gca_900096595.Microbotrium_Intermedium_Assembly.dna.toplevel.fa \
     -I MI949596sorted_dedup_RGI-splitreads.bam \
     -O raw_varients949596.vcf \
     -bamout bamout.bam
```

## 13. Select variant types (SNPs / indels)

```bash
gatk SelectVariants \
     -R Microbotryum_intermedium_gca_900096595.Microbotrium_Intermedium_Assembly.dna.toplevel.fa \
     -V raw_varients949596.vcf \
     --select-type-to-include SNP \
     -O raw_snps.vcf

gatk SelectVariants \
     -R Microbotryum_intermedium_gca_900096595.Microbotrium_Intermedium_Assembly.dna.toplevel.fa \
     -V raw_varients949596.vcf \
     --select-type-to-include INDEL \
     -O raw_indels.vcf
```

## 14. Filter variants (hard filtering)

SNPs:

```bash
gatk VariantFiltration \
     -V raw_snps979899.vcf \
     -filter "QD<2.0"              --filter-name "QD2" \
     -filter "QUAL<30.0"           --filter-name "QUAL30" \
     -filter "SOR>3.0"             --filter-name "SOR3" \
     -filter "FS>60.0"             --filter-name "FS60" \
     -filter "MQ<40.0"             --filter-name "MQ40" \
     -filter "MQRankSum<-12.5"     --filter-name "MQRankSum-12.5" \
     -filter "ReadPosRankSum<-8.0" --filter-name "ReadPosRankSum-8" \
     -O raw_snps-filtered.vcf
```

Indels:

```bash
gatk VariantFiltration \
     -V raw_indels979899.vcf \
     -filter "QD<2.0"               --filter-name "QD2" \
     -filter "QUAL<30.0"            --filter-name "QUAL30" \
     -filter "FS>200.0"             --filter-name "FS200" \
     -filter "ReadPosRankSum<-20.0" --filter-name "ReadPosRankSum-20" \
     -O raw_indels979899-filtered.vcf
```

## 15. Build the snpEff annotation database

```bash
java -jar snpEff.jar build -gtf22 -v genome_m_intermedium
```

## 16. Annotate variants

```bash
java -Xmx8g -jar snpEff.jar genome_m_intermedium \
     raw_snps-filtered.vcf > raw_snps-filtered.ann.vcf

java -Xmx8g -jar snpEff.jar genome_m_intermedium \
     raw_indels-filtered.vcf > raw_indels-filtered.ann.vcf
```

## 17. Convert VCF to BED and annotate against gene coordinates

Uses [BEDOPS](https://bedops.readthedocs.io/).

**Step 1** — convert the annotation GTF (`genes.gtf`) file into a BED file:

```bash
gtf2bed < genes.gtf > sorted-genes.gtf.bed
```

**Step 2** — convert the final filtered/annotated VCF into a BED file:

```bash
vcf2bed < raw_snps979899-filtered.ann.vcf | sort-bed - > sites979899.bed
```

**Step 3** — combine the two BED files from the previous steps into a final annotated
BED file:

```bash
bedmap --echo --echo-map-id --delim '\t' \
       sites979899.bed sorted-genes.gtf.bed > annotatedSites979899.bed
```

Output: `annotatedSites979899.bed`, used to generate the figures in
[04_visualization_R](../04_visualization_R/plots.R).
