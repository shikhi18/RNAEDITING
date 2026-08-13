# Step 1: Quality Trimming & Adapter Removal (fastp)

This step takes raw paired-end reads, processes them for quality improvement, and
outputs cleaned paired-end reads in `.fastq.gz` format using
[fastp](https://github.com/OpenGene/fastp).

## Basic run

```bash
fastp -i 18134R-216-01_S94_L001_R1_001.fastq.gz \
      -I 18134R-216-01_S94_L001_R2_001.fastq.gz \
      -o R1_001.fastq.gz \
      -O R2_001.fastq.gz
```

## Set a minimum read length

Sets the minimum length for reads to be retained. Any read shorter than 36 bases after
trimming is discarded (`-l`).

```bash
fastp -i 18134R-216-01_S94_L001_R1_001.fastq.gz \
      -I 18134R-216-01_S94_L001_R2_001.fastq.gz \
      -o R1_001.fastq.gz \
      -O R2_001.fastq.gz \
      -l 36
```

## Auto-detect adapters for paired-end data

Enables automatic detection of adapter sequences for paired-end sequencing data. fastp
identifies and trims adapter sequences based on the overlap between paired reads
(`--detect_adapter_for_pe`).

```bash
fastp -i 18134R-216-01_S94_L001_R1_001.fastq.gz \
      -I 18134R-216-01_S94_L001_R2_001.fastq.gz \
      -o R1_001.fastq.gz \
      -O R2_001.fastq.gz \
      -l 36 \
      --detect_adapter_for_pe
```

## Per-base quality trimming from both ends

Useful for trimming low-quality sequences that tend to occur at the ends of reads
(`--cut_right`, `--cut_front`).

```bash
fastp -i 18134R-216-01_S94_L001_R1_001.fastq.gz \
      -I 18134R-216-01_S94_L001_R2_001.fastq.gz \
      -o R1_001.fastq.gz \
      -O R2_001.fastq.gz \
      -l 36 \
      --cut_right --cut_front
```

## Final command (with global quality cutting enabled)

In fastp, `--cut_right`, `--cut_front`, and `-c` (base correction for overlapping
paired-end reads) are used together to fine-tune the trimming process based on quality.

```bash
fastp -i 18134R-216-01_S94_L001_R1_001.fastq.gz \
      -I 18134R-216-01_S94_L001_R2_001.fastq.gz \
      -o R1_001.fastq.gz \
      -O R2_001.fastq.gz \
      -l 36 \
      --cut_right --cut_front -c
```

Next step: [02_alignment_STAR](../02_alignment_STAR/README.md)
