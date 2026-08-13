# RNA Editing Analysis — *Microbotryum intermedium*

Scripts and pipeline documentation for an RNA-editing study in *Microbotryum
intermedium* (anther-smut fungus), covering RNA-seq processing, GATK-based
variant calling, and the R scripts used to generate the paper/thesis figures.

Sample groups referenced in the pipeline docs: single-sample runs (e.g.
sample `94`) and merged-replicate runs (`94`/`95`/`96` and `97`/`98`/`99`).
Figure scripts reference multiple developmental/life-cycle stages
(`p1A1`, `p1A2`, `6P`, `6D`, `MIA1`, `MIA2`, and a `Mitate`/infection stage).

## Repository contents

This repo is a flat collection of pipeline-step write-ups and the R scripts
used for figures — it is **not** organized into numbered subfolders, even
though some of the docs below cross-reference folder names like
`02_alignment_STAR/`. Those references are historical and don't correspond
to real paths in this repo; use the file list below instead.

| File | Purpose |
|---|---|
| [`Quality Trimming & Adapter Removal (fastp).md`](Quality%20Trimming%20%26%20Adapter%20Removal%20%28fastp%29.md) | Step 1 — read trimming/QC with `fastp` |
| [`fastp codes.odt`](fastp%20codes.odt) | Raw terminal transcript of the actual `fastp` + STAR indexing/alignment commands run, with real output stats |
| [`Genome Indexing & Alignment (STAR).md`](Genome%20Indexing%20%26%20Alignment%20%28STAR%29.md) | Step 2 — GFF3→GTF conversion, STAR genome indexing, and alignment |
| [`GATK.md`](GATK.md) | Step 3 — duplicate marking, `SplitNCigarReads`, `HaplotypeCaller`, hard filtering, `snpEff` annotation, and BED conversion via BEDOPS |
| [`Fig9ACODE(stackplot).R`](Fig9ACODE%28stackplot%29.R) | Figure 9A — stacked bar plot of editing-site functional categories by stage |
| [`FIG 9B-RCODE.R`](FIG%209B-RCODE.R) | Figure 9B — bar plot of amino-acid substitution frequency by species/haploid group |
| [`FIG10B-RCODE.R`](FIG10B-RCODE.R) | Figure 10B — bar plot of `log2FoldChange` by gene function |
| [`MVLG_MILATE_UNIQUE GENES_stage_wordcloud.R`](MVLG_MILATE_UNIQUE%20GENES_stage_wordcloud.R) | Word cloud of MVLG genes uniquely edited in the Mitate stage, by functional category |
| [`MvSup Common Genes Word Cloud Plot.R`](MvSup%20Common%20Genes%20Word%20Cloud%20Plot.R) | Word cloud of MvSup genes common to all four stages, by functional category |
| [`hap1-RCODE.R`](hap1-RCODE.R) | Haplotype 1 (`fig1ef`) bar plot of amino-acid substitution frequency, faceted by haploid group |

## Pipeline overview

1. **Trim & QC reads** — `fastp`, see the fastp doc above. Minimum length 36 bp,
   automatic paired-end adapter detection, quality-based cutting from both
   ends, and overlap-based base correction (`-l 36 --detect_adapter_for_pe
   --cut_right --cut_front -c`).
2. **Align to genome** — STAR, see the STAR doc above. Genome annotation is
   converted from GFF3 to GTF with `gffread` before indexing.
3. **Call and filter variants** — GATK best-practices RNA-seq short-variant
   workflow (`MarkDuplicates` → `AddOrReplaceReadGroups` →
   `SplitNCigarReads` → `HaplotypeCaller` → `SelectVariants` →
   `VariantFiltration` → `snpEff` → BEDOPS `vcf2bed`/`bedmap`), producing an
   `annotatedSites*.bed` file per sample group. See `GATK.md`.
4. **Summarize sites into figures** — the R scripts read pre-summarized
   Excel workbooks (`NEW2.xlsx`, `fig1efmay16th.xlsx`) with per-gene,
   per-stage editing frequencies and functional annotations, and produce the
   bar plots and word clouds. **These workbooks are not included in this
   repo**, and the step that turns `annotatedSites*.bed` into those
   per-gene/per-stage summary tables (site-level frequency/coverage
   filtering, functional annotation curation) is not scripted here — see
   Known limitations below.

## Tools used

`fastp` (v0.20.1), `STAR` (v2.7.10a), `cufflinks`/`gffread`, `samtools`,
`GATK` (4.3.0.0), `picard`, `snpEff`, `BEDOPS`, and R packages `ggplot2`,
`readxl`, `openxlsx`, `reshape2`, `wordcloud`, `RColorBrewer`.

## Known limitations

- `GATK.md`'s merged-replicate steps 13–17 have inconsistent intermediate
  file names across steps (e.g. step 13 writes `raw_snps.vcf`, but step 14
  reads `raw_snps979899.vcf`) — likely commands from the two replicate
  groups (`94`/`95`/`96` vs. `97`/`98`/`99`) were interleaved when this was
  written up. Treat the commands as a reference for parameters/thresholds
  rather than a literally copy-pasteable script until the file names are
  reconciled.
- The doc-to-doc links originally pointed at a numbered-folder layout
  (`02_alignment_STAR/`, `03_variant_calling_GATK/`, `04_visualization_R/`)
  that doesn't exist in this repo; `04_visualization_R/plots.R` in
  particular is referenced from `GATK.md` but was never added.
- There is no script or documentation in this repo for the step between
  `annotatedSites*.bed` and the summary Excel workbooks the R scripts
  consume — i.e., the actual RNA-editing site filtering/curation logic
  (e.g. selecting A-to-G/T-to-C conversions, coverage or frequency
  thresholds, functional-annotation curation) isn't captured in code here.
