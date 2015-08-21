###########################################################
# You may be able to get this to work with shell          #
# I tend to lean bash and have used that environment here #
###########################################################
SHELL := /bin/bash

#####################################
# Key macros used extensively below #
# and values loaded from files      #
#####################################
RELEASE := 79
PY_CPU_FILE := num_ngg2_cpus
MULTI_CPUS := $(shell cat ${PY_CPU_FILE})
SPECIES_FILE := species
SPECIES := $(shell cat ${SPECIES_FILE})

#################################################################
# Try to keep directory clean with subdirs for specific outputs #
#################################################################
DIRS := grna_out gc_content logs gtfs fasta R_output figures benchmark
GTFFILES := $(addprefix gtfs/, $(addsuffix .gtf, $(SPECIES)))
FASTAFILES := $(addprefix fasta/, $(addsuffix .fa, $(SPECIES)))
GCFILES := $(addprefix gc_content/, $(addsuffix _gcPerc.txt, $(SPECIES)))

######################################################################
# Used in summaries and figures (exhaustive search, all gRNA starts) #
######################################################################
GRNA_FULL_OUTPUT := $(addprefix grna_out/, $(addsuffix _full_grnas.csv, $(SPECIES)))

###################################
# Benchmark only, canonical sites #
###################################
BENCHMARK_BLOCK_1CPU := $(addprefix benchmark/, $(addsuffix _benchmark_block_1cpu_grnas.log, $(SPECIES)))
BENCHMARK_FULL_1CPU := $(addprefix benchmark/, $(addsuffix _benchmark_full_1cpu_grnas.log, $(SPECIES)))
BENCHMARK_BLOCK_MULTICPU := $(addprefix benchmark/, $(addsuffix _benchmark_block_multicpu_grnas.log, $(SPECIES)))
BENCHMARK_FULL_MULTICPU := $(addprefix benchmark/, $(addsuffix _benchmark_full_multicpu_grnas.log, $(SPECIES)))

###################################################
# Basic analysis output i.e. gene / gRNA overlaps #
###################################################
RFILEWITHPATH = $(addprefix R_output/, all_species_counts.csv all_species_pam.csv all_species_strand.csv all_species_gstart_pam.csv all_species_gstart_strand.csv all_species_all_gRNA_biotypes.csv all_species_unique_gRNA_biotypes.csv all_species_gstart_gRNA_biotypes.csv all_species_gstart_unique_gRNA_biotypes.csv)

all: $(DIRS) $(FASTAFILES) $(GCFILES) $(GTFFILES) $(GRNA_FULL_OUTPUT) $(RFILEWITHPATH) $(BENCHMARK_BLOCK_1CPU) $(BENCHMARK_FULL_1CPU) $(BENCHMARK_BLOCK_MULTICPU) $(BENCHMARK_FULL_MULTICPU)
	@echo Macros used in this makefile
	@echo DIRS: $(DIRS)
	@echo FASTA: $(FASTAFILES)
	@echo GC: $(GCFILES)
	@echo GTF: $(GTFFILES)
	@echo RFILES: $(RFILEWITHPATH)
	@echo SPECIES: $(SPECIES)
	@echo NGG2 multicpu: $(MULTI_CPUS)

########################################
# Makes sure correct directories exist #
########################################
$(DIRS):
	mkdir -p $(DIRS)

#######################################
# Calculate GC content of each genome #
#######################################
gc_content/%_gcPerc.txt: fasta/%.fa
	./fasta_gc_content_calc.py $< --output_file $@ 1>logs/$*_gc.log 2>logs/$*_gc.log
	
gtfs/%.gtf:
	wget --quiet ftp://ftp.ensembl.org/pub/release-${RELEASE}/gtf/$*/*.gtf.gz -O gtfs/$*.gtf.gz
	gunzip gtfs/$*.gtf.gz
	
grna_out/%_full_grnas.csv: fasta/%.fa
	/usr/bin/time -f "%P CPU\n%M Max mem (kb)\n%S CPU-seconds\n%e real-time" ./ngg2_regex_unique_multproc.py --allowNoncanonical --cores ${MULTI_CPUS} $< --outputFile $@ 1>logs/$*_full_multicpu.log 2>logs/$*_full_multicpu_error.log
	
$(RFILEWITHPATH): $(GRNA_FULL_OUTPUT) $(GTFFILES) $(FASTAFILES) ngg2_analysis.R
	R --vanilla < ./ngg2_analysis.R 1>logs/R_annotations.log 2>logs/R_annotations.log

benchmark/%_benchmark_full_multicpu_grnas.log: fasta/%.fa
	/usr/bin/time -f "%e real-time" ./ngg2_regex_unique_multproc.py --cores ${MULTI_CPUS} $< --outputFile tmp.grna 1>tmp.log 2>$@
	rm tmp.grna tmp.log

benchmark/%_benchmark_full_1cpu_grnas.log:  fasta/%.fa
	/usr/bin/time -f "%e real-time" ./ngg2_regex_unique_multproc.py --cores 1 $< --outputFile tmp.grna 1>tmp.log 2>$@
	rm tmp.grna tmp.log

benchmark/%_benchmark_block_multicpu_grnas.log:  fasta/%.fa
	/usr/bin/time -f "%e real-time" ./ngg2_regex_unique_multproc.py --blockScan --cores ${MULTI_CPUS} $< --outputFile tmp.grna 1>tmp.log 2>$@
	rm tmp.grna tmp.log
	
benchmark/%_benchmark_block_1cpu_grnas.log:  fasta/%.fa
	/usr/bin/time -f "%e real-time" ./ngg2_regex_unique_multproc.py --blockScan --cores 1 $< --outputFile tmp.grna 1>tmp.log 2>$@
	rm tmp.grna tmp.log
	
fasta/%.fa:
	wget --quiet ftp://ftp.ensembl.org/pub/release-${RELEASE}/fasta/$*/dna/*dna_sm.primary_assembly.fa.gz -O fasta/$*.fa.gz
	
	if [ ! -s fasta/$*.fa.gz ]; \
	then wget --quiet ftp://ftp.ensembl.org/pub/release-${RELEASE}/fasta/$*/dna/*dna_sm.toplevel.fa.gz -O fasta/$*.fa.gz; \
	fi;
	gunzip fasta/$*.fa.gz
	