# viridi_taxonomy_workflow

I developed this Snakemake workflow to replace the ad hoc scripts I used to generate the carbohydrate-active enzyme (CAZyme), KEGG, and taxonomic annotation files for my 2020 master's thesis. It brings those reads-based metagenomic analyses into a single codebase with explicit inputs, outputs, and reference paths.

The workflow includes DIAMOND searches against CAZyme, KEGG, and RefSeq databases, along with MEGAN rules for taxonomic assignment and per-sample count export. The root Snakefile defaults to Bowtie2 alignment of paired-end reads. The taxonomy Snakefile defaults to exporting per-sample taxonomic counts from precomputed DIAMOND archives.
