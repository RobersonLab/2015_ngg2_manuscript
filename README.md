# 2015_ngg2_manuscript
Repository for scripts associated with manuscript.

## Directory layout
The top level directory contains the scripts necessary to repeat the analysis.

The figures directory contains the images generated during the final analysis using R Markdown.

## Running the analysis

### Requirements
* R, packages will prompt for install if needed
* make
* wget
* samtools
* R Studio with knitr
* python, along with the [ngg2 program](https://github.com/RobersonLab/ngg2)
* the scripts from this repository

### Downloading required genomes
The makefile handles downloading the genomes, creating the FASTA index for each genome, downloading the gene annotations, calculating genome GC content for each organism, and running ngg2 on each genome. It should be run first.

```bash
make -f ngg2_datageneration.make 1>prep_data.log 2>&1
```

### Annotating gRNA overlaps with genes and gather meta information
The annotation of gRNAs and genes uses the Bioconductor GenomicRanges package. The data preparation should download all files necessary to annotate the data. Annotation is performed in R using ngg2_analysis.R.

The script is set by default to not allow parallel execution. However, parallel execution will greatly decrease runtime. Adjust the PARALLEL variable to TRUE and set the number of CPUs to run the ldply and ddply sections in parallel. doMC is registered as the backend, but that could be altered to suit your preferences.

```bash
R --vanilla < ngg2_analysis.R 1>ngg.log 2>&1
```

### Generate summary information, tables, and figures with R Markdown
The summarizing R Markdown file ngg2_plots_and_summaries.RMd should be in the same directory as the analysis output. To recreate the analysis, simple open the Markdown file in R Studio and knit the output to your preferred file format.
