# Step 2: Genome Indexing & Alignment (STAR)

## Convert genome annotation from GFF3 to GTF

STAR requires a GTF annotation file. This converts the GFF3 annotation to GTF format
(using the [cufflinks](http://cole-trapnell-lab.github.io/cufflinks/) package's
`gffread`).

```bash
gffread Microbotryum_intermedium_gca_900096595.Microbotrium_Intermedium_Assembly.56.gff3 \
        -T \
        -o cufflinks_M.intermedium_output.gtf
```

## Generate the genome index

This command generates a set of genome index files in the `ref/` directory. These index
files are used later to align RNA-seq reads to the reference genome using STAR in
alignment mode.

```bash
STAR --runMode genomeGenerate \
     --genomeDir ref/ \
     --genomeFastaFiles Microbotryum_intermedium_gca_900096595.Microbotrium_Intermedium_Assembly.dna.toplevel.fa \
     --sjdbGTFfile cufflinks_M.intermedium_output.gtf
```

## Align reads

```bash
STAR --runMode alignReads \
     --genomeDir ref/ \
     --outSAMtype BAM SortedByCoordinate \
     --readFilesIn R1_02_S95_001.fastq R2_02_S95_001.fastq
```

Output: a coordinate-sorted BAM file (e.g. `Aligned.sortedByCoord.out.bam`), used as
input for [Step 3: variant calling with GATK](../03_variant_calling_GATK/README.md).
