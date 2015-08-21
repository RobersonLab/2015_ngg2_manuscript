# 2015_ngg2_manuscript
Repository for scripts associated with manuscript.

## Directory layout
The top level directory contains the scripts necessary to repeat the analysis.

The figures directory contains the images generated during the final analysis using R Markdown.

## Running the analysis

### Requirements
* Linux + bash
* R, packages should prompt for install if needed
* make
* wget
* R Studio with knitr (to knit the RMd file)
* python, along with the [ngg2 program](https://github.com/RobersonLab/ngg2)
* the scripts from this repository

### Downloading required genomes
The makefile handles the entire analysis from file downloads to annoation of sites. There are several parameters you may wish to edit.

* species - this file contains a list of all species to download from Ensembl. Add or remove as you wish.
* num_ngg2_cpus - sets how many processors to run ngg2
* num_rscript_cpus - sets how many processors to run R analysis (may be memory intensive)
* ngg2_analysis.R - you may wish to change the registered parallel backend for R to use (make sure to deactivate at end of R script)

### Running analysis
```bash
nohup make -f ngg2_datageneration.make 1>ngg2_paper.log 2>&1 &
```

## Associated paper data
There are two versions of the data associated with this paper.

- Version 1
> This was run with the early versions of ngg2 and ran a **block** scan. That means that if two gRNA sites overlap one another on the same strand, only the first is reported. Only 6 model species were tested.

- Version 2
> This version was run with a more efficient, multiprocessing ngg2 version using **exhaustive** scanning to identify all 3'GG gRNA sites. Two additional species, African elephan and pigs, were added in this analysis.

| Data | Type | Preprint |
| :--- | :--- | :---- |
| [Version 1 fileset](http://dx.doi.org/10.6084/m9.figshare.1371077) | Block scan | [ngg2 v1](https://dx.doi.org/10.7287/peerj.preprints.969v1) |
| [Version 2 fileset](http://dx.doi.org/10.6084/m9.figshare.1515944) | Exhaustive scan | |
