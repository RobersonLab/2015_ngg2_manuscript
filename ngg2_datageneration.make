RELEASE = 79
SPECIES = homo_sapiens mus_musculus danio_rerio drosophila_melanogaster caenorhabditis_elegans saccharomyces_cerevisiae
FASTAFILES = $(addsuffix .fa, $(SPECIES))
FASTAINDEXES = $(addsuffix .fai, $(FASTAFILES))
GRNAFILES = $(addsuffix _grnas.csv, $(SPECIES))
GTFFILES = $(addsuffix .gtf, $(SPECIES))
GCFILES = $(addsuffix _gcPerc.txt, $(SPECIES))

all: $(FASTAFILES) $(FASTAINDEXES) $(GCFILES) $(GTFFILES) $(GRNAFILES)

%_gcPerc.txt: %.fa
	./gcCalc.py $*.fa --outputFile $*_gcPerc.txt > $*_gc.log

%.fa.fai: %.fa
	samtools faidx $*.fa

%_grnas.csv:  %.fa %.fa.fai
	ngg2.py --allContigs $*.fa --outputFile $*_grnas.csv > $*_grna.log

%.gtf:
	wget --quiet ftp://ftp.ensembl.org/pub/release-${RELEASE}/gtf/$*/*.gtf.gz -O $*.gtf.gz
	gunzip $*.gtf.gz

homo_sapiens.fa:
	wget --quiet ftp://ftp.ensembl.org/pub/release-${RELEASE}/fasta/homo_sapiens/dna/*dna_sm.primary_assembly.fa.gz -O homo_sapiens.fa.gz
	gunzip homo_sapiens.fa.gz

mus_musculus.fa:
	wget --quiet ftp://ftp.ensembl.org/pub/release-${RELEASE}/fasta/mus_musculus/dna/*dna_sm.primary_assembly.fa.gz -O mus_musculus.fa.gz
	gunzip mus_musculus.fa.gz

danio_rerio.fa:
	wget --quiet ftp://ftp.ensembl.org/pub/release-${RELEASE}/fasta/danio_rerio/dna/*dna_sm.toplevel.fa.gz -O danio_rerio.fa.gz
	gunzip danio_rerio.fa.gz

drosophila_melanogaster.fa:
	wget --quiet ftp://ftp.ensembl.org/pub/release-${RELEASE}/fasta/drosophila_melanogaster/dna/*dna_sm.toplevel.fa.gz -O drosophila_melanogaster.fa.gz
	gunzip drosophila_melanogaster.fa.gz

caenorhabditis_elegans.fa:
	wget --quiet ftp://ftp.ensembl.org/pub/release-${RELEASE}/fasta/caenorhabditis_elegans/dna/*dna_sm.toplevel.fa.gz -O caenorhabditis_elegans.fa.gz
	gunzip caenorhabditis_elegans.fa.gz

saccharomyces_cerevisiae.fa:
	wget --quiet ftp://ftp.ensembl.org/pub/release-${RELEASE}/fasta/saccharomyces_cerevisiae/dna/*dna_sm.toplevel.fa.gz -O saccharomyces_cerevisiae.fa.gz
	gunzip saccharomyces_cerevisiae.fa.gz

